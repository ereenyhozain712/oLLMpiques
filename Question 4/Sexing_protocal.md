# High-Throughput Sexing of *Cannabis sativa* Using REM16 Expression

## 1. Marker Choice
REM16 separates the sexes completely: females ~10.0–11.1, males ~7.95–9.4 (p < 1e-16), with a gap near 9.7. FT1 shows no difference (p = 0.69), so it is excluded.

**REM16 alone is used as the classifier.**

## 2. Sample Collection
1. Take a ~4 mm leaf punch from each seedling.
2. Place it in a barcoded 96-well plate with RNA-stabilizing buffer (or on dry ice).
3. Record plant ID and plate position in a sample sheet.

## 3. Sample Processing
1. Homogenize tissue with a bead mill.
2. Extract RNA using a magnetic-bead, plate-based kit on a liquid-handling robot.
3. Run one-step RT-qPCR in 384-well plates, duplexing REM16 with a stable reference gene.
4. Include controls on every plate: **known male**, **known female**, and **no-template**.

## 4. Sex Determination
- Calculate **ΔCt = Ct(reference) − Ct(REM16)**.
- Females (higher REM16) have a lower ΔCt than males.
- Set the threshold at the midpoint of the gap between groups, calibrated on controls.
- Flag samples near the threshold, or with failed reference amplification, for retesting.

## 5. Scaling Up
- Automated pipetting and 384-well plates allow **1,000+ plants per day**.
- An R/Python script reads qPCR output, applies the threshold, and returns a sex-call table.
- Testing at the seedling stage lets growers remove unwanted males before flowering, saving space and cost.

## 6. Validation
- Confirm calls against phenotype or an XY-linked DNA marker in a pilot batch, and test across tissues, stages, and cultivars.
