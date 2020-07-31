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
