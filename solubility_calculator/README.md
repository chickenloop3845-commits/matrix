# Methylphenidate Analog Water Solubility Calculator

A comprehensive Python tool for calculating aqueous solubility of methylphenidate and its analogs using multiple parameter-based methods.

## Features

### Calculation Methods

1. **General Solubility Equation (GSE)** - Yalkowsky method
   - `log S = 0.5 - 0.01(MP - 25) - log P`
   - Uses melting point and partition coefficient

2. **Modified GSE with Structural Corrections**
   - Adds corrections for H-bonding, aromaticity, flexibility
   - Improved accuracy for drug-like molecules

3. **Abraham Solvation Parameters (LFER)**
   - Uses E, S, A, B, V descriptors
   - Linear free energy relationship approach

4. **Hansen Solubility Parameters**
   - 3D solubility parameter space (δD, δP, δH)
   - Distance-based solubility prediction

5. **Group Contribution Method**
   - Fragment-based calculation
   - Additive contributions from molecular groups

6. **QSPR Model**
   - Multiple linear regression
   - Uses logP, MW, TPSA, HBD, HBA, rotatable bonds

7. **pH-Dependent Solubility**
   - Henderson-Hasselbalch equation
   - Accounts for ionization state

8. **Salt Form Predictions**
   - HCl, sulfate, phosphate, tartrate salts
   - Empirical salt factors

9. **Temperature Dependence**
   - van't Hoff equation
   - Solubility at multiple temperatures

10. **Consensus Prediction**
    - Average of multiple methods
    - Statistical confidence assessment

## Compounds Included

All 20 methylphenidate analogs:

1. Dexmethylphenidate (d-threo-MPH)
2. L-threo-Methylphenidate
3. Ethylphenidate
4. 3,4-Dichloroethylphenidate
5. 4-Fluoromethylphenidate (4F-MPH)
6. 4-Fluoroethylphenidate
7. 3,4-Dichloromethylphenidate (3,4-CTMP)
8. Isopropylphenidate (IPH)
9. Propylphenidate
10. 4-Methylmethylphenidate
11. HDMP-28
12. HDEP-28
13. Serdexmethylphenidate
14. 4-Chloromethylphenidate
15. 3-Methylmethylphenidate
16. Cyclopropylphenidate
17. N-Ethyl-methylphenidate
18. tert-Butylphenidate
19. 4-Bromomethylphenidate
20. 2-Fluoromethylphenidate

## Usage

### Basic Usage

```python
from methylphenidate_solubility import (
    METHYLPHENIDATE_ANALOGS,
    SolubilityCalculator,
    SaltForm
)

# Get a compound
compound = METHYLPHENIDATE_ANALOGS["dexmethylphenidate"]

# Create calculator
calc = SolubilityCalculator(compound)

# Calculate using GSE method
result = calc.calculate_gse()
print(f"Solubility: {result.solubility_mg_mL:.4f} mg/mL")
```

### pH-Dependent Solubility

```python
# Calculate at different pH values
for pH in [1.0, 4.0, 7.0, 7.4]:
    result = calc.calculate_ph_dependent(pH)
    print(f"pH {pH}: {result.solubility_mg_mL:.4f} mg/mL")
```

### Salt Form Comparison

```python
# Compare different salt forms
for salt in SaltForm:
    result = calc.calculate_salt_solubility(salt)
    print(f"{salt.value}: {result.solubility_mg_mL:.2f} mg/mL")
```

### Consensus Prediction

```python
# Get consensus from all methods
consensus = calc.calculate_consensus()
print(f"Mean: {consensus['consensus'].solubility_mg_mL:.4f} mg/mL")
print(f"Std Dev: {consensus['statistics']['std_log_S']:.3f} log units")
```

### Compare Multiple Compounds

```python
from methylphenidate_solubility import compare_compounds

compounds = ["dexmethylphenidate", "ethylphenidate", "4_fluoromethylphenidate"]
print(compare_compounds(compounds))
```

### Generate Full Report

```python
from methylphenidate_solubility import generate_full_report

report = generate_full_report("dexmethylphenidate")
print(report)
```

## Running the Calculator

```bash
cd solubility_calculator
python methylphenidate_solubility.py
```

## Output Format

Each calculation returns a `SolubilityResult` object containing:

- `compound_name`: Name of the compound
- `method`: Calculation method used
- `solubility_mol_L`: Solubility in mol/L
- `solubility_mg_mL`: Solubility in mg/mL
- `log_solubility`: log10 of solubility (mol/L)
- `temperature_C`: Temperature in Celsius
- `pH`: pH value (if applicable)
- `salt_form`: Salt form used
- `confidence`: Prediction confidence (high/medium/low)
- `notes`: Additional calculation details

## Molecular Descriptors

Each compound includes:

- Molecular weight
- LogP (octanol-water partition coefficient)
- Melting point
- pKa
- Hydrogen bond donors/acceptors
- Topological polar surface area (TPSA)
- Rotatable bonds
- Aromatic rings
- Molar volume
- Molar refractivity
- Abraham parameters (E, S, A, B, V)
- Hansen parameters (δD, δP, δH)

## Theoretical Background

### General Solubility Equation (GSE)

The GSE relates aqueous solubility to two key properties:
- **Melting point**: Higher MP = lower solubility (crystal lattice energy)
- **LogP**: Higher lipophilicity = lower aqueous solubility

### Abraham Parameters

Linear free energy relationship using:
- **E**: Excess molar refraction (polarizability)
- **S**: Dipolarity/polarizability
- **A**: Hydrogen bond acidity (donor strength)
- **B**: Hydrogen bond basicity (acceptor strength)
- **V**: McGowan characteristic volume

### Hansen Solubility Parameters

Three-dimensional approach:
- **δD**: Dispersion forces
- **δP**: Polar forces
- **δH**: Hydrogen bonding

Solubility decreases with increasing "distance" from water in Hansen space.

### pH-Dependent Solubility

For basic compounds like methylphenidate:
```
S_total = S_0 × (1 + 10^(pKa - pH))
```

At low pH, the compound is ionized and more soluble.

## References

1. Yalkowsky, S.H. & Valvani, S.C. (1980) J. Pharm. Sci. 69, 912-922.
2. Abraham, M.H. (1993) Chem. Soc. Rev. 22, 73-83.
3. Hansen, C.M. (2007) Hansen Solubility Parameters: A User's Handbook.
4. Lipinski, C.A. et al. (2001) Adv. Drug Deliv. Rev. 46, 3-26.

## License

For research and educational purposes only.

## Disclaimer

This calculator provides theoretical predictions based on established models. Experimental validation is recommended for any practical applications. The solubility values are estimates and may differ from actual experimental measurements.
