*ENTRNA: A Framework to Predict RNA Foldability*
=========================================================================

The RNA foldability problem is defined as: 
> Given a pair of RNA sequence and secondary structure, what is the likelihood that the given RNA sequence would fold into the given secondary structure. 

RNA foldability prediction fundamentally differs from RNA secondary structure prediction and the RNA inverse folding problem, as the later ones only require RNA sequences or secondary structure as a single input. RNA foldability prediction will require both sequence and secondary structure to be provided. However, RNA foldability could be applied to solve RNA secondary structure prediction and the RNA inverse folding problem, since both of them require accurate accessment for the fitness of a pair of RNA sequence and secodnary structure.

ENTRNA is a data-driven framework for RNA foldability prediction. In addition, it is able to evaluate both pseudoknot-free and pseudoknotted RNAs.

## Requirements

* Python 2.7
* Vienna RNA python2 bindings 
* numpy (== 1.13.3)
* pandas (== 0.19.2)
* sklearn (== 0.19.1)
* scipy ( == 0.19.0)

## Usage

ENTRNA takes sequence and secondary structure as inputs. 
To deal with both pseudoknot-free and pseudoknotted RNAs, it uses base-pairing array to represent the RNA secondary structure. 
For example, a pseudoknot-free secondary structure in dot-bracket notation `.(((...))).` could be represented as `[0,10,9,8,0,0,0,3,2,1,0]`.

Two components are provided in `ENTRNA`: `ENTRNA_train` and `ENTRNA_prediction`.

### ENTRNA_train

`ENTRNA_train` is for training ENTRNA framework. 

#### Example
- script:
```shell
python ENTRNA_train.py --is_pseudoknot_free y --real_rna_path ./util/RNASTRAND_pseudoknot_free_feature.csv --simulation_rna_path ./util/RNASTRAND_extract_feature_pseudoknot_free/
```

- output:
```
Training ENTRNA for pseudoknot-free RNAs
ENTRNA training accuracy: 0.832517140059
```



### ENTRNA_predict
`ENTRNA_predict` is for RNA foldability prediction. It is pretrained based on RNAs extracted from [RNASTRAND database](http://www.rnasoft.ca/strand/)

#### Example:
- script:
```shell
python ENTRNA_predict.py --seq_file pseudoknotted_seq.txt --str_file pseudoknotted_str.txt 
```

- output:
```
RNA sequence:
GGCGCGGCACCGUCCGCGGAACAAACGG
RNA secondary structure:
[0, 0, 18, 17, 16, 15, 14, 0, 0, 28, 27, 26, 0, 7, 6, 5, 4, 3, 0, 0, 0, 0, 0, 0, 0, 12, 11, 10]
This is pseudoknotted RNA
Foldability: 0.947906134401
```


