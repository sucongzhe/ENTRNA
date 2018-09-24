ENTRNA: *A Framework to Predict RNA Foldability*
=========================================================================

The RNA foldability problem is defined as: 
> Given a pair of RNA sequence and secondary structure, what is the likelihood that the given RNA sequence would fold into the given secondary structure. 

RNA foldability prediction fundamentally differs from RNA secondary structure prediction and the RNA inverse folding problem, as the later ones only require RNA sequences or secondary structure as a single input. RNA foldability prediction will require both sequence and secondary structure to be provided. However, RNA foldability could be applied to solve RNA secondary structure prediction and the RNA inverse folding problem, since both of them require accurate accessment for the fitness of a pair of RNA sequence and secodnary structure.

ENTRNA is a data-driven framework for RNA foldability prediction. In addition, it is able to evaluate both pseudoknot-free and pseudoknotted RNAs.

Requirements
------------

* Python 2.7
* Vienna RNA python2 bindings 
* numpy (== 1.13.3)
* pandas (== 0.19.2)
* sklearn (== 0.19.1)
* scipy ( == 0.19.0)

Usage
-----

ENTRNA takes sequence and secondary structure as inputs. To deal with both pseudoknot-free and pseudoknotted RNAs, it uses base-pairing array to represent the RNA secondary structure. For example, a pseudoknot-free secondary structure in dot-bracket notation `.(((...))).` could be represented as `[0,10,9,8,0,0,0,3,2,1,0]`.

`ENTRNA` is the main program of this package. It extracts features for the pair of RNA sequence and secondary structure, and trains a classification model, and predicts the RNA foldability. 

Example:
```python
from ENTRNA import entrna_main
seq = 'GGACUCAGUAAUAUGCUUUGGAAACGAAGCUUACAAAAUGGAGUCC'
bp = [46, 45, 44, 43, 42, 41, 0, 0, 0, 0, 0, 0, 0, 0, 30, 29, 28, 27, 26, 25, 0, 0, 0, 0, 20, 19, 18, 17, 16, 15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 5, 4, 3, 2, 1]
foldability = entrna_main(seq,bp)
```
