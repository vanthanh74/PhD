## Usage

```
python3 main.py [-h] [--requested-rank REQUESTED_RANKS [REQUESTED_RANKS ...]]
               [--cores NBS_CORES [NBS_CORES ...]] [--single-rank]
               [--disp-cols]
               matrix [matrix ...]

positional arguments:
  matrix                matrix path

optional arguments:
  -h, --help            show this help message and exit
  --requested-rank REQUESTED_RANKS [REQUESTED_RANKS ...], -k REQUESTED_RANKS [REQUESTED_RANKS ...]
  --cores NBS_CORES [NBS_CORES ...], -c NBS_CORES [NBS_CORES ...]
                        n for square splitting, nxm for rectangular splitting,
                        append 'c' at then end for column then row gathering
                        instead of nested squares gathering
  --single-rank         single_rank
  --disp-cols           Display selected columns at each step
```

# Example

```bash
python3 main.py --single-rank ../../code/matrix_generation/csv/phillips.txt -c 1x2,2x1 -k 10
```
python3 main.py --single-rank csvMatrices/baart.txt -c 1x2,2x1 -k 10
