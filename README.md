ENTRNA: A Framework to Predict RNA Foldability
=========================================================================

Requirements
------------

* Python 2.7
* Vienna RNA python bindings 
* numpy (== 1.13.3)
* pandas (== 0.19.2)
* sklearn (== 0.19.1)
* scipy ( == 0.19.0)

Usage
-----

ENTRNA takes sequence and base-pairing array as inputs. 
Both pseudoknot-free and pseudoknotted RNAs could be evaluated.


Here is the sample code:
```python
import ENTRNA
seq = 'GGACUCAGUAAUAUGCUUUGGAAACGAAGCUUACAAAAUGGAGUCC'
bp = [46, 45, 44, 43, 42, 41, 0, 0, 0, 0, 0, 0, 0, 0, 30, 29, 28, 27, 26, 25, 0, 0, 0, 0, 20, 19, 18, 17, 16, 15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 5, 4, 3, 2, 1]
foldability = ENTRNA.entrna_main(seq,bp)
```
