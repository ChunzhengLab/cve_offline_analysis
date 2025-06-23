# Feeddown Dispose Script

This script applies feeddown correction to Proton data using Lambda secondary data.

## Formula

The correction uses the following formula:
```
O_lp = (O_data - (1 - fp)O_lpsec)/fp
```

Where:
- `O_data`: Original Proton data (delta or gamma)
- `O_lpsec`: Lambda secondary data (from Lambda files)
- `fp`: Primary fraction (hardcoded average of p and ap primary fractions)
- `O_lp`: Corrected Proton data

## Error Propagation

Error is calculated using:
```
σ_Olp = (1 / f_p) * sqrt[σ_D² + (1 - f_p)² * σ_S² + (O_lpsec - O_lp)² * σ_f²]
```

## Usage

### Process all tasks
```bash
python feeddown_dispose.py
```

### Process specific task
```bash
python feeddown_dispose.py --task default
python feeddown_dispose.py -t ChiMax2
```

### Custom directories
```bash
python feeddown_dispose.py --proton-dir /path/to/proton --lambda-dir /path/to/lambda --output-dir /path/to/output
```

## Input Data

### Proton Data
- Location: `../dataset_merger/Merged/Proton/finalise_[task].csv`
- Format: `centrality,diff_type,diff_bin,pair_type,delta,delta_err,gamma,gamma_err`

### Lambda Data
- Location: `../dataset_merger/Merged/Lambda/finalise_[task].csv`
- Fallback: `../dataset_merger/Merged/Lambda/finalise_default.csv` if task-specific file not found
- Format: `centrality,pair_type,delta,delta_err,gamma,gamma_err[,resolution]`

## Primary Fraction Values

Hardcoded values used for calculation:
- p primary total: `0.8489101154204324 ± 0.12451480395843792`
- ap primary total: `0.854837091993265 ± 0.1336355938419346`
- fp = average of p and ap primary values

## Processing Logic

1. Load Proton data for each task
2. Load corresponding Lambda data (or default if not found)
3. For each row with `pair_type` ∈ {SS, OS}:
   - Apply feeddown correction to both delta and gamma
   - Use Lambda data with matching centrality and pair_type
4. Recalculate `Del` rows as `OS - SS` using corrected values
5. Save corrected data to output directory

## Output

Corrected files are saved to:
`./Proton/finalise_feeddown_dispose_[task].csv`

## Command Line Options

- `--task`, `-t`: Specify task name to process
- `--proton-dir`: Custom Proton data directory
- `--lambda-dir`: Custom Lambda data directory  
- `--output-dir`: Custom output directory
- `--help`, `-h`: Show help message
- `--version`: Show version information

## Examples

```bash
# Process all tasks
python feeddown_dispose.py

# Process only default task
python feeddown_dispose.py -t default

# Use custom directories
python feeddown_dispose.py --proton-dir ./custom/proton --output-dir ./results

# Get help
python feeddown_dispose.py --help
```