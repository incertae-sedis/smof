2.22.3 [2023-03-03]

 * Fix bug in `subseq`

2.22.2 [2022-07-18]

 * Fix bug in `consensus`
 * Fix bug in `subseq --color`

2.22.1 [2021-07-19]

  * fix instant death bug in `cut`
  * add tests for `cut`

2.22.0 [2021-05-14]

  * add --delimiter argument to clean for removing descriptions 

2.21.1 [2021-02-15]

  * fix bug in translation of gapped DNA sequences

2.21.0 [2020-09-14]

  * modularize functions for import into other Python packages
  * make the subseq annotation optional ("|subseq(i..j)")

2.20.0 [2020-09-03]

  * draft modularization

2.19.0 [2020-07-29]

  * change `smof uniq --final-header` to `--first-header` and keep first header
    rather than last
  * add `smof uniq --removed=FILE` option that writes removed entries to FILE

2.18.0 [2020-03-02]

  * allow lowercase 'X' to be replace with 'n' in clean standard function
  * raise warning when character type is not specified in clean standard function
  * add `--cds` option to `translate` that returns DNA coding sequences

2.17.0 [2019-12-31]

  * force order to be retained in uniq -f

  * add 'standardize' option to clean that replaces 'X' in DNA with 'N' and
    standardizes gaps from `[._-]` to just '-' 
