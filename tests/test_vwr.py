import io
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
from synth import write_vwr

from eeg2edf import edfcommon
from eeg2edf.micromed import convert as vwr2edf
from eeg2edf.micromed import vwr

try:
    import edfio
except ImportError:
    edfio = None


class VwrReaderTest(unittest.TestCase):
    def test_header_order_calibration_and_notes(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "record.vwr"
            values = write_vwr(path)
            header = vwr.read_header(path)

            self.assertEqual(header.start.isoformat(), "2024-01-02T03:04:05")
            self.assertEqual([channel.name for channel in header.channels], ["B-REF", "A-REF"])
            self.assertEqual([(n.frame, n.description) for n in header.notes],
                             [(0, "start"), (2, "onset")])
            frames = next(vwr.iter_frames(header, header.n_samples))
            np.testing.assert_array_equal(frames, values)
            np.testing.assert_allclose(vwr.calibrated(frames, header.channels),
                                       [[1, 1], [2, 2], [3, 3], [4, 4]])

    def test_all_libvwr_sample_widths(self):
        with tempfile.TemporaryDirectory() as directory:
            for width in (1, 2, 4):
                path = Path(directory) / f"record-{width}.vwr"
                values = write_vwr(path, width)
                header = vwr.read_header(path)
                np.testing.assert_array_equal(next(vwr.iter_frames(header, 4)), values)

    def test_rejects_incomplete_frames(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "broken.vwr"
            write_vwr(path)
            with open(path, "ab") as fh:
                fh.write(b"\0")
            with self.assertRaisesRegex(ValueError, "trailing"):
                vwr.read_header(path)


class VwrMontageTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.path = Path(self.directory.name) / "record.vwr"
        write_vwr(self.path)
        self.header = vwr.read_header(self.path)

    def tearDown(self):
        self.directory.cleanup()

    def test_montages_are_found_by_name_not_segment_position(self):
        # MONTAGE sits between ORDER and LABCOD in the synthetic file.
        self.assertEqual([m.name for m in self.header.montages], ["Bipolar", "Ground"])
        self.assertEqual(self.header.montage("Bipolar").notch, 1)

    def test_input_zero_resolves_to_the_channels_own_ground(self):
        traces = self.header.montage("Bipolar").traces
        self.assertEqual([t.label for t in traces], ["A-B", "B-REF"])
        self.assertEqual((traces[1].active, traces[1].reference, traces[1].reference_index),
                         ("B", "REF", 0))

    def test_history_names_the_recorded_montage(self):
        self.assertEqual(self.header.recorded_montage, "Ground")

    def test_montage_pairs_map_onto_order_positions(self):
        # ORDER is (LABCOD 2, LABCOD 1), so B is position 0 and A is position 1.
        pairs = vwr2edf.montage_pairs(self.header, self.header.montage("Bipolar"))
        self.assertEqual(pairs, [("A-B", 1, 0), ("B-REF", 0, None)])

    def test_montage_output_labels_and_values(self):
        output = Path(self.directory.name) / "montage.edf"
        montage = self.header.montage("Bipolar")
        self.assertEqual(vwr2edf.convert(self.header, output, "X X X X", True, montage), 2)
        meta = json.loads(output.with_suffix(".json").read_text())
        self.assertEqual([c["edf_label"] for c in meta["channels"]], ["A-B", "B-REF"])
        self.assertEqual(meta["montage_applied"], "Bipolar")
        self.assertEqual([c["derived"] for c in meta["channels"]], [True, False])
        if edfio:
            edf = edfio.read_edf(io.BytesIO(output.read_bytes()))
            self.assertEqual(edf.labels, ("A-B", "B-REF"))
            # Both channels calibrate to 1..4, so their difference is flat zero.
            # A-B spans twice a 16-bit channel's range in int16, so it is exact
            # only to one output LSB -- the sidecar records how big that is.
            lsb = meta["channels"][0]["resolution_uv_per_lsb"]
            np.testing.assert_allclose(edf.signals[0].data, [0, 0, 0, 0], atol=lsb)
            np.testing.assert_allclose(edf.signals[1].data, [1, 2, 3, 4], atol=0.01)

    def test_unknown_montage_is_rejected(self):
        self.assertIsNone(self.header.montage("nope"))


class VwrMarkerTest(unittest.TestCase):
    def test_every_marker_area_is_read_and_sorted(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "record.vwr"
            write_vwr(path)
            header = vwr.read_header(path)
            self.assertEqual(
                [(m.frame, m.label, m.source, m.end_frame) for m in header.markers],
                [(0, "Flag 1", "FLAGS", 2),
                 (0, "start", "NOTE", None),
                 (1, "Seizure", "EVENT A", 3),
                 (1, "Trigger 7", "TRIGGER", None),
                 (2, "onset", "NOTE", None)],
            )
            self.assertEqual([(p.frame, p.original_frame) for p in header.parts],
                             [(0, 10), (2, 20)])

    def test_intervals_become_annotations_with_a_duration(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "record.vwr"
            write_vwr(path)
            header = vwr.read_header(path)
            events = vwr2edf._events(header)
            self.assertEqual(events[0], (0.0, "Flag 1", 0.5))
            self.assertEqual(events[2], (0.25, "Seizure", 0.5))
            self.assertIn(b"+0\x150.5\x14Flag 1", edfcommon._tal(0.0, "Flag 1", 0.5))


class VwrEdfTest(unittest.TestCase):
    def test_edf_and_sidecar(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "record.vwr"
            output = Path(directory) / "record.edf"
            write_vwr(source)
            header = vwr.read_header(source)
            vwr2edf.convert(header, output, "X X X X", sidecar=True)

            raw = output.read_bytes()
            self.assertEqual(raw[0:8].rstrip(), b"0")
            self.assertIn(b"onset", raw)
            self.assertIn(b"Trigger 7", raw)
            header_size = 256 * (header.order + 2)
            samples = np.frombuffer(raw, dtype="<i2",
                                    count=header.order * header.frequency, offset=header_size)
            self.assertEqual(samples.size, 8)
            self.assertTrue(np.all(samples >= -32768))

            metadata = json.loads(output.with_suffix(".json").read_text())
            self.assertEqual(metadata["schema"], edfcommon.SCHEMA)
            self.assertEqual(metadata["clip"]["sfreq_hz"], 4)
            self.assertEqual(metadata["source"]["format"], "micromed-vwr")
            self.assertIsNone(metadata["montage_applied"])
            self.assertEqual(metadata["events"][0]["label"], "Flag 1")
            self.assertEqual(metadata["events"][0]["duration_s"], 0.5)
            self.assertEqual(len(metadata["segments"]), 2)
            self.assertNotIn("patient_name", metadata)
            self.assertEqual(json.dumps(metadata["patient"]),
                             json.dumps({"sex": None, "dob": None, "age_at_recording": None}))
            if edfio:
                edf = edfio.read_edf(io.BytesIO(output.read_bytes()))
                self.assertEqual(edf.labels, ("B-REF", "A-REF"))
                np.testing.assert_allclose(edf.signals[0].data, [1, 2, 3, 4], atol=1.1)
                np.testing.assert_allclose(edf.signals[1].data, [1, 2, 3, 4], atol=1.1)
                self.assertEqual(len(edf.annotations), 5)

    def test_last_record_is_padded_at_physical_zero(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "record.vwr"
            output = Path(directory) / "record.edf"
            write_vwr(source)
            header = vwr.read_header(source)
            # 4 samples at 4 Hz is exactly one record; ask for a longer record
            # so the tail is padded.
            object.__setattr__(header, "frequency", 8)
            vwr2edf.convert(header, output, "X X X X", sidecar=False)
            if edfio:
                edf = edfio.read_edf(io.BytesIO(output.read_bytes()))
                np.testing.assert_allclose(edf.signals[0].data[4:], 0, atol=0.11)


if __name__ == "__main__":
    unittest.main()
