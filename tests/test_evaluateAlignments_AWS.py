import logging
import unittest
import time
import evaluateAlignments_AWS as evalScript

class LoggerTestCase(unittest.TestCase):
    def test_check_timeout(self):
        def timeout_afterTime(sec):
            t_end = time.time() + sec
            while time.time() < t_end:
                2 * 2
        self.assertEqual(logging.debug('not killed yet!'), evalScript.check_timeout(timeout_afterTime(5), timeout=3))


if __name__ == '__main__':
    unittest.main()
[{},
 {'prob': '100.0', 'eval': '5.8E-91', 'pval': '4.2E-95', 'hhscore': '665.8', 'aligned_cols': '330', 'q_range': '1-330',
  't_range': '1-330', 'match md5': 'f853053b9c86b1f72f32db6bc3d945b0', 'pdbCode': '12as_B', 'pdbRange': '4-330',
  'pdbMatchLength': '327', 'cathCodes': ['3.30.930.10'], 'identities': '100', 'similarity': '1.480',
  'maxcluster': {'avrg': {}, 'max': {}, 'min': {}, 'range': {}}}]
