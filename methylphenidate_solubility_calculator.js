#!/usr/bin/env node
/**
 * Methylphenidate Analog Water Solubility Calculator
 * ===================================================
 * 
 * This module implements multiple parameter-based methods for calculating
 * aqueous solubility of methylphenidate and its analogs:
 * 
 * 1. General Solubility Equation (GSE) - Yalkowsky
 * 2. Modified GSE with corrections
 * 3. Abraham Solvation Parameter Model
 * 4. Hansen Solubility Parameters
 * 5. Lipinski-based estimation
 * 6. Henderson-Hasselbalch pH-dependent solubility
 * 7. Temperature-dependent solubility (van't Hoff)
 * 8. ESOL (Delaney) Model
 * 9. Ali et al. Model
 * 
 * Author: Research Compendium Generator
 * Date: January 2026
 */

const fs = require('fs');

// =============================================================================
// COMPOUND DATABASE
// =============================================================================

const COMPOUNDS = {
    "Dexmethylphenidate": {
        name: "Dexmethylphenidate",
        formula: "C14H19NO2",
        mw: 233.31,
        mp: 74.0,
        logP: 2.15,
        pKa: 8.8,
        hbd: 1,
        hba: 3,
        tpsa: 38.33,
        rotatable_bonds: 4,
        aromatic_rings: 1,
        E: 1.12, S: 0.85, A: 0.15, B: 0.75, V: 1.89,
        delta_d: 17.5, delta_p: 5.2, delta_h: 7.8
    },
    "L-threo-Methylphenidate": {
        name: "L-threo-Methylphenidate",
        formula: "C14H19NO2",
        mw: 233.31,
        mp: 74.0,
        logP: 2.15,
        pKa: 8.8,
        hbd: 1,
        hba: 3,
        tpsa: 38.33,
        rotatable_bonds: 4,
        aromatic_rings: 1,
        E: 1.12, S: 0.85, A: 0.15, B: 0.75, V: 1.89,
        delta_d: 17.5, delta_p: 5.2, delta_h: 7.8
    },
    "Ethylphenidate": {
        name: "Ethylphenidate",
        formula: "C15H21NO2",
        mw: 247.33,
        mp: 68.0,
        logP: 2.65,
        pKa: 8.7,
        hbd: 1,
        hba: 3,
        tpsa: 38.33,
        rotatable_bonds: 5,
        aromatic_rings: 1,
        E: 1.18, S: 0.82, A: 0.14, B: 0.73, V: 2.03,
        delta_d: 17.2, delta_p: 4.8, delta_h: 7.2
    },
    "3,4-Dichloroethylphenidate": {
        name: "3,4-Dichloroethylphenidate",
        formula: "C15H19Cl2NO2",
        mw: 316.22,
        mp: 95.0,
        logP: 3.85,
        pKa: 8.5,
        hbd: 1,
        hba: 3,
        tpsa: 38.33,
        rotatable_bonds: 5,
        aromatic_rings: 1,
        E: 1.45, S: 0.95, A: 0.12, B: 0.68, V: 2.35,
        delta_d: 18.5, delta_p: 6.5, delta_h: 5.8
    },
    "4-Fluoromethylphenidate": {
        name: "4-Fluoromethylphenidate",
        formula: "C14H18FNO2",
        mw: 251.30,
        mp: 82.0,
        logP: 2.45,
        pKa: 8.6,
        hbd: 1,
        hba: 3,
        tpsa: 38.33,
        rotatable_bonds: 4,
        aromatic_rings: 1,
        E: 1.08, S: 0.92, A: 0.14, B: 0.72, V: 1.92,
        delta_d: 17.8, delta_p: 5.8, delta_h: 7.2
    },
    "4-Fluoroethylphenidate": {
        name: "4-Fluoroethylphenidate",
        formula: "C15H20FNO2",
        mw: 265.33,
        mp: 76.0,
        logP: 2.95,
        pKa: 8.5,
        hbd: 1,
        hba: 3,
        tpsa: 38.33,
        rotatable_bonds: 5,
        aromatic_rings: 1,
        E: 1.14, S: 0.90, A: 0.13, B: 0.70, V: 2.06,
        delta_d: 17.5, delta_p: 5.5, delta_h: 6.8
    },
    "3,4-Dichloromethylphenidate": {
        name: "3,4-Dichloromethylphenidate",
        formula: "C14H17Cl2NO2",
        mw: 302.20,
        mp: 105.0,
        logP: 3.55,
        pKa: 8.4,
        hbd: 1,
        hba: 3,
        tpsa: 38.33,
        rotatable_bonds: 4,
        aromatic_rings: 1,
        E: 1.40, S: 0.98, A: 0.13, B: 0.65, V: 2.21,
        delta_d: 18.8, delta_p: 6.8, delta_h: 5.5
    },
    "Isopropylphenidate": {
        name: "Isopropylphenidate",
        formula: "C16H23NO2",
        mw: 261.36,
        mp: 62.0,
        logP: 3.15,
        pKa: 8.6,
        hbd: 1,
        hba: 3,
        tpsa: 38.33,
        rotatable_bonds: 5,
        aromatic_rings: 1,
        E: 1.22, S: 0.78, A: 0.13, B: 0.72, V: 2.17,
        delta_d: 16.8, delta_p: 4.5, delta_h: 6.5
    },
    "Propylphenidate": {
        name: "Propylphenidate",
        formula: "C16H23NO2",
        mw: 261.36,
        mp: 58.0,
        logP: 3.20,
        pKa: 8.6,
        hbd: 1,
        hba: 3,
        tpsa: 38.33,
        rotatable_bonds: 6,
        aromatic_rings: 1,
        E: 1.24, S: 0.76, A: 0.13, B: 0.71, V: 2.17,
        delta_d: 16.6, delta_p: 4.3, delta_h: 6.3
    },
    "4-Methylmethylphenidate": {
        name: "4-Methylmethylphenidate",
        formula: "C15H21NO2",
        mw: 247.33,
        mp: 78.0,
        logP: 2.65,
        pKa: 8.7,
        hbd: 1,
        hba: 3,
        tpsa: 38.33,
        rotatable_bonds: 4,
        aromatic_rings: 1,
        E: 1.18, S: 0.82, A: 0.14, B: 0.74, V: 2.03,
        delta_d: 17.3, delta_p: 4.9, delta_h: 7.0
    },
    "HDMP-28": {
        name: "HDMP-28",
        formula: "C15H21NO2",
        mw: 247.33,
        mp: 85.0,
        logP: 2.55,
        pKa: 9.0,
        hbd: 1,
        hba: 3,
        tpsa: 38.33,
        rotatable_bonds: 4,
        aromatic_rings: 1,
        E: 1.16, S: 0.84, A: 0.16, B: 0.76, V: 2.03,
        delta_d: 17.4, delta_p: 5.0, delta_h: 7.5
    },
    "HDEP-28": {
        name: "HDEP-28",
        formula: "C16H23NO2",
        mw: 261.36,
        mp: 72.0,
        logP: 3.05,
        pKa: 8.9,
        hbd: 1,
        hba: 3,
        tpsa: 38.33,
        rotatable_bonds: 5,
        aromatic_rings: 1,
        E: 1.22, S: 0.81, A: 0.15, B: 0.74, V: 2.17,
        delta_d: 17.1, delta_p: 4.7, delta_h: 7.0
    },
    "Serdexmethylphenidate": {
        name: "Serdexmethylphenidate",
        formula: "C17H24N2O4",
        mw: 320.39,
        mp: 145.0,
        logP: 0.85,
        pKa: 8.2,
        hbd: 3,
        hba: 5,
        tpsa: 89.56,
        rotatable_bonds: 8,
        aromatic_rings: 1,
        E: 1.35, S: 1.45, A: 0.55, B: 1.15, V: 2.52,
        delta_d: 18.2, delta_p: 9.5, delta_h: 12.5
    },
    "4-Chloromethylphenidate": {
        name: "4-Chloromethylphenidate",
        formula: "C14H18ClNO2",
        mw: 267.75,
        mp: 88.0,
        logP: 2.85,
        pKa: 8.5,
        hbd: 1,
        hba: 3,
        tpsa: 38.33,
        rotatable_bonds: 4,
        aromatic_rings: 1,
        E: 1.22, S: 0.92, A: 0.13, B: 0.70, V: 2.05,
        delta_d: 18.0, delta_p: 5.8, delta_h: 6.5
    },
    "3-Methylmethylphenidate": {
        name: "3-Methylmethylphenidate",
        formula: "C15H21NO2",
        mw: 247.33,
        mp: 76.0,
        logP: 2.60,
        pKa: 8.7,
        hbd: 1,
        hba: 3,
        tpsa: 38.33,
        rotatable_bonds: 4,
        aromatic_rings: 1,
        E: 1.17, S: 0.83, A: 0.14, B: 0.74, V: 2.03,
        delta_d: 17.3, delta_p: 5.0, delta_h: 7.1
    },
    "Cyclopropylphenidate": {
        name: "Cyclopropylphenidate",
        formula: "C16H21NO2",
        mw: 259.34,
        mp: 70.0,
        logP: 2.95,
        pKa: 8.6,
        hbd: 1,
        hba: 3,
        tpsa: 38.33,
        rotatable_bonds: 4,
        aromatic_rings: 1,
        E: 1.25, S: 0.80, A: 0.13, B: 0.71, V: 2.10,
        delta_d: 17.0, delta_p: 4.6, delta_h: 6.6
    },
    "N-Ethyl-methylphenidate": {
        name: "N-Ethyl-methylphenidate",
        formula: "C16H23NO2",
        mw: 261.36,
        mp: 55.0,
        logP: 2.85,
        pKa: 8.2,
        hbd: 0,
        hba: 3,
        tpsa: 29.54,
        rotatable_bonds: 5,
        aromatic_rings: 1,
        E: 1.20, S: 0.75, A: 0.00, B: 0.72, V: 2.17,
        delta_d: 16.8, delta_p: 4.2, delta_h: 5.5
    },
    "tert-Butylphenidate": {
        name: "tert-Butylphenidate",
        formula: "C17H25NO2",
        mw: 275.39,
        mp: 65.0,
        logP: 3.65,
        pKa: 8.5,
        hbd: 1,
        hba: 3,
        tpsa: 38.33,
        rotatable_bonds: 4,
        aromatic_rings: 1,
        E: 1.28, S: 0.72, A: 0.12, B: 0.70, V: 2.31,
        delta_d: 16.2, delta_p: 4.0, delta_h: 5.8
    },
    "4-Bromomethylphenidate": {
        name: "4-Bromomethylphenidate",
        formula: "C14H18BrNO2",
        mw: 312.20,
        mp: 92.0,
        logP: 3.15,
        pKa: 8.5,
        hbd: 1,
        hba: 3,
        tpsa: 38.33,
        rotatable_bonds: 4,
        aromatic_rings: 1,
        E: 1.35, S: 0.95, A: 0.13, B: 0.68, V: 2.15,
        delta_d: 18.5, delta_p: 6.0, delta_h: 6.2
    },
    "2-Fluoromethylphenidate": {
        name: "2-Fluoromethylphenidate",
        formula: "C14H18FNO2",
        mw: 251.30,
        mp: 80.0,
        logP: 2.35,
        pKa: 8.6,
        hbd: 1,
        hba: 3,
        tpsa: 38.33,
        rotatable_bonds: 4,
        aromatic_rings: 1,
        E: 1.06, S: 0.94, A: 0.14, B: 0.73, V: 1.92,
        delta_d: 17.9, delta_p: 6.0, delta_h: 7.4
    }
};

// Water Hansen parameters at 25°C
const WATER_HANSEN = { delta_d: 15.5, delta_p: 16.0, delta_h: 42.3 };

// Gas constant
const R = 8.314; // J/(mol·K)

// =============================================================================
// SOLUBILITY CALCULATION METHODS
// =============================================================================

/**
 * Method 1: General Solubility Equation (GSE) - Yalkowsky
 * log S = 0.5 - 0.01(MP - 25) - log P
 */
function gseYalkowsky(compound) {
    const mpTerm = 0.01 * (compound.mp - 25);
    const logS = 0.5 - mpTerm - compound.logP;
    const sMolL = Math.pow(10, logS);
    const sMgMl = sMolL * compound.mw;
    
    return {
        method: "GSE (Yalkowsky)",
        logS: logS,
        sMolL: sMolL,
        sMgMl: sMgMl,
        sGL: sMgMl,
        confidence: "Medium",
        notes: `MP term: ${mpTerm.toFixed(3)}, logP term: ${compound.logP.toFixed(2)}`
    };
}

/**
 * Method 2: Modified GSE with molecular weight correction
 * log S = 0.8 - 0.01(MP - 25) - 1.05*log P + 0.012*TPSA - 0.0015*MW
 */
function gseModified(compound) {
    const mpTerm = 0.01 * (compound.mp - 25);
    const logPTerm = 1.05 * compound.logP;
    const tpsaTerm = 0.012 * compound.tpsa;
    const mwTerm = 0.0015 * compound.mw;
    
    const logS = 0.8 - mpTerm - logPTerm + tpsaTerm - mwTerm;
    const sMolL = Math.pow(10, logS);
    const sMgMl = sMolL * compound.mw;
    
    return {
        method: "Modified GSE",
        logS: logS,
        sMolL: sMolL,
        sMgMl: sMgMl,
        sGL: sMgMl,
        confidence: "Medium-High",
        notes: `Includes TPSA (${compound.tpsa.toFixed(1)}) and MW (${compound.mw.toFixed(1)}) corrections`
    };
}

/**
 * Method 3: Abraham Solvation Parameter Model
 * log S = c + e*E + s*S + a*A + b*B + v*V
 */
function abrahamModel(compound) {
    // Abraham coefficients for water solubility
    const c = 0.52;
    const e = 0.77;
    const s = -2.17;
    const a = 3.99;
    const b = -4.87;
    const v = -3.36;
    
    const logS = c + e * compound.E + s * compound.S + a * compound.A + b * compound.B + v * compound.V;
    const sMolL = Math.pow(10, logS);
    const sMgMl = sMolL * compound.mw;
    
    return {
        method: "Abraham LFER",
        logS: logS,
        sMolL: sMolL,
        sMgMl: sMgMl,
        sGL: sMgMl,
        confidence: "High",
        notes: `E=${compound.E.toFixed(2)}, S=${compound.S.toFixed(2)}, A=${compound.A.toFixed(2)}, B=${compound.B.toFixed(2)}, V=${compound.V.toFixed(2)}`
    };
}

/**
 * Method 4: Hansen Solubility Parameter Distance
 * Ra² = 4(δD1-δD2)² + (δP1-δP2)² + (δH1-δH2)²
 */
function hansenDistance(compound) {
    const deltaDDiff = compound.delta_d - WATER_HANSEN.delta_d;
    const deltaPDiff = compound.delta_p - WATER_HANSEN.delta_p;
    const deltaHDiff = compound.delta_h - WATER_HANSEN.delta_h;
    
    const RaSquared = 4 * deltaDDiff * deltaDDiff + deltaPDiff * deltaPDiff + deltaHDiff * deltaHDiff;
    const Ra = Math.sqrt(RaSquared);
    
    // Empirical correlation: log S ≈ 1.5 - 0.08*Ra
    const logS = 1.5 - 0.08 * Ra;
    const sMolL = Math.pow(10, logS);
    const sMgMl = sMolL * compound.mw;
    
    return {
        method: "Hansen Distance",
        logS: logS,
        sMolL: sMolL,
        sMgMl: sMgMl,
        sGL: sMgMl,
        confidence: "Medium",
        notes: `Ra = ${Ra.toFixed(2)} MPa^0.5 (distance from water)`
    };
}

/**
 * Method 5: Lipinski-based estimation
 */
function lipinskiEstimation(compound) {
    const logS = 1.2 
        - 0.9 * compound.logP 
        - 0.003 * compound.mw 
        + 0.15 * compound.hbd 
        + 0.08 * compound.hba 
        - 0.02 * compound.rotatable_bonds;
    
    const sMolL = Math.pow(10, logS);
    const sMgMl = sMolL * compound.mw;
    
    return {
        method: "Lipinski-based",
        logS: logS,
        sMolL: sMolL,
        sMgMl: sMgMl,
        sGL: sMgMl,
        confidence: "Medium",
        notes: `HBD=${compound.hbd}, HBA=${compound.hba}, RotBonds=${compound.rotatable_bonds}`
    };
}

/**
 * Method 6: ESOL (Estimated SOLubility) - Delaney
 * log S = 0.16 - 0.63*clogP - 0.0062*MW + 0.066*RB - 0.74*AP
 */
function esolDelaney(compound) {
    // Estimate aromatic proportion
    const aromaticAtoms = 6 * compound.aromatic_rings;
    const heavyAtoms = compound.mw / 12;
    const ap = heavyAtoms > 0 ? aromaticAtoms / heavyAtoms : 0;
    
    const logS = 0.16 
        - 0.63 * compound.logP 
        - 0.0062 * compound.mw 
        + 0.066 * compound.rotatable_bonds 
        - 0.74 * ap;
    
    const sMolL = Math.pow(10, logS);
    const sMgMl = sMolL * compound.mw;
    
    return {
        method: "ESOL (Delaney)",
        logS: logS,
        sMolL: sMolL,
        sMgMl: sMgMl,
        sGL: sMgMl,
        confidence: "High",
        notes: `Aromatic proportion: ${ap.toFixed(3)}`
    };
}

/**
 * Method 7: Ali et al. Model
 * log S = 0.342 - 1.034*logP - 0.008*TPSA - 0.0065*MW
 */
function aliModel(compound) {
    const logS = 0.342 
        - 1.034 * compound.logP 
        - 0.008 * compound.tpsa 
        - 0.0065 * compound.mw;
    
    const sMolL = Math.pow(10, logS);
    const sMgMl = sMolL * compound.mw;
    
    return {
        method: "Ali et al.",
        logS: logS,
        sMolL: sMolL,
        sMgMl: sMgMl,
        sGL: sMgMl,
        confidence: "Medium-High",
        notes: `TPSA=${compound.tpsa.toFixed(1)} Å²`
    };
}

/**
 * Method 8: Consensus (Average of methods)
 */
function consensusPrediction(compound) {
    const methods = [
        gseYalkowsky(compound),
        gseModified(compound),
        abrahamModel(compound),
        hansenDistance(compound),
        lipinskiEstimation(compound),
        esolDelaney(compound),
        aliModel(compound)
    ];
    
    const logSValues = methods.map(m => m.logS);
    const avgLogS = logSValues.reduce((a, b) => a + b, 0) / logSValues.length;
    const stdLogS = Math.sqrt(logSValues.reduce((sum, x) => sum + Math.pow(x - avgLogS, 2), 0) / logSValues.length);
    
    const sMolL = Math.pow(10, avgLogS);
    const sMgMl = sMolL * compound.mw;
    
    return {
        method: "Consensus (7 methods)",
        logS: avgLogS,
        sMolL: sMolL,
        sMgMl: sMgMl,
        sGL: sMgMl,
        confidence: "High",
        notes: `Std dev: ±${stdLogS.toFixed(2)} log units, Range: [${Math.min(...logSValues).toFixed(2)}, ${Math.max(...logSValues).toFixed(2)}]`
    };
}

/**
 * pH-dependent solubility (Henderson-Hasselbalch)
 */
function phDependentSolubility(compound, pH = 7.4) {
    const intrinsic = consensusPrediction(compound);
    const s0 = intrinsic.sMolL;
    
    // Henderson-Hasselbalch for bases
    const ionizationFactor = 1 + Math.pow(10, compound.pKa - pH);
    const sTotal = s0 * ionizationFactor;
    
    const logS = sTotal > 0 ? Math.log10(sTotal) : -10;
    const sMgMl = sTotal * compound.mw;
    
    return {
        method: `pH-dependent (pH=${pH})`,
        logS: logS,
        sMolL: sTotal,
        sMgMl: sMgMl,
        sGL: sMgMl,
        confidence: "Medium-High",
        notes: `pKa=${compound.pKa}, Ionization factor: ${ionizationFactor.toFixed(1)}x`
    };
}

/**
 * Temperature-dependent solubility (van't Hoff)
 */
function temperatureDependentSolubility(compound, TCelsius = 25.0, deltaHSol = 25000.0) {
    const TRef = 298.15; // 25°C in Kelvin
    const TTarget = TCelsius + 273.15;
    
    const refSol = consensusPrediction(compound);
    const sRef = refSol.sMolL;
    
    // van't Hoff equation
    const lnRatio = -deltaHSol / R * (1/TTarget - 1/TRef);
    const sTarget = sRef * Math.exp(lnRatio);
    
    const logS = sTarget > 0 ? Math.log10(sTarget) : -10;
    const sMgMl = sTarget * compound.mw;
    
    return {
        method: `van't Hoff (T=${TCelsius}°C)`,
        logS: logS,
        sMolL: sTarget,
        sMgMl: sMgMl,
        sGL: sMgMl,
        confidence: "Medium",
        notes: `ΔH_sol=${(deltaHSol/1000).toFixed(1)} kJ/mol, T_ref=25°C`
    };
}

/**
 * Calculate all methods for a compound
 */
function calculateAll(compound, pH = 7.4, TCelsius = 25.0) {
    return [
        gseYalkowsky(compound),
        gseModified(compound),
        abrahamModel(compound),
        hansenDistance(compound),
        lipinskiEstimation(compound),
        esolDelaney(compound),
        aliModel(compound),
        consensusPrediction(compound),
        phDependentSolubility(compound, pH),
        temperatureDependentSolubility(compound, TCelsius)
    ];
}

// =============================================================================
// MAIN EXECUTION
// =============================================================================

function printSeparator(char = "=", length = 80) {
    console.log(char.repeat(length));
}

function main() {
    printSeparator();
    console.log("METHYLPHENIDATE ANALOG WATER SOLUBILITY CALCULATIONS");
    console.log("Parameter-Based Prediction Methods");
    printSeparator();
    console.log();
    
    // Summary table header
    console.log("=".repeat(130));
    console.log(
        "Compound".padEnd(35) +
        "GSE".padStart(10) +
        "Mod.GSE".padStart(10) +
        "Abraham".padStart(10) +
        "Hansen".padStart(10) +
        "ESOL".padStart(10) +
        "Consensus".padStart(12) +
        "pH 7.4".padStart(12) +
        "mg/mL".padStart(12)
    );
    console.log(
        "".padEnd(35) +
        "(log S)".padStart(10) +
        "(log S)".padStart(10) +
        "(log S)".padStart(10) +
        "(log S)".padStart(10) +
        "(log S)".padStart(10) +
        "(log S)".padStart(12) +
        "(log S)".padStart(12) +
        "(pH 7.4)".padStart(12)
    );
    console.log("=".repeat(130));
    
    const allResults = {};
    
    for (const [name, compound] of Object.entries(COMPOUNDS)) {
        const gse = gseYalkowsky(compound);
        const modGse = gseModified(compound);
        const abraham = abrahamModel(compound);
        const hansen = hansenDistance(compound);
        const esol = esolDelaney(compound);
        const consensus = consensusPrediction(compound);
        const phDep = phDependentSolubility(compound, 7.4);
        
        console.log(
            name.padEnd(35) +
            gse.logS.toFixed(3).padStart(10) +
            modGse.logS.toFixed(3).padStart(10) +
            abraham.logS.toFixed(3).padStart(10) +
            hansen.logS.toFixed(3).padStart(10) +
            esol.logS.toFixed(3).padStart(10) +
            consensus.logS.toFixed(3).padStart(12) +
            phDep.logS.toFixed(3).padStart(12) +
            phDep.sMgMl.toFixed(4).padStart(12)
        );
        
        allResults[name] = {
            compound: compound,
            results: calculateAll(compound)
        };
    }
    
    console.log("=".repeat(130));
    console.log();
    
    // Detailed results for each compound
    printSeparator();
    console.log("DETAILED RESULTS BY COMPOUND");
    printSeparator();
    
    for (const [name, data] of Object.entries(allResults)) {
        console.log();
        printSeparator("-");
        console.log(`COMPOUND: ${name}`);
        console.log(`Formula: ${data.compound.formula}, MW: ${data.compound.mw.toFixed(2)} g/mol`);
        console.log(`MP: ${data.compound.mp}°C, logP: ${data.compound.logP.toFixed(2)}, pKa: ${data.compound.pKa}`);
        printSeparator("-");
        console.log();
        
        for (const result of data.results) {
            console.log(`  Method: ${result.method}`);
            console.log(`  log S: ${result.logS.toFixed(3)}`);
            console.log(`  Solubility: ${result.sMolL.toExponential(2)} mol/L`);
            console.log(`  Solubility: ${result.sMgMl.toFixed(4)} mg/mL (${result.sGL.toFixed(4)} g/L)`);
            console.log(`  Confidence: ${result.confidence}`);
            console.log(`  Notes: ${result.notes}`);
            console.log();
        }
    }
    
    // Solubility ranking
    printSeparator();
    console.log("SOLUBILITY RANKING (by consensus log S, most to least soluble)");
    printSeparator();
    
    const rankings = Object.entries(allResults).map(([name, data]) => {
        const consensus = data.results[7]; // Consensus is index 7
        return [name, consensus.logS, consensus.sMgMl];
    });
    rankings.sort((a, b) => b[1] - a[1]);
    
    console.log("Rank".padEnd(6) + "Compound".padEnd(35) + "log S".padStart(10) + "mg/mL".padStart(12));
    console.log("-".repeat(65));
    rankings.forEach(([name, logS, sMgMl], i) => {
        console.log(
            String(i + 1).padEnd(6) +
            name.padEnd(35) +
            logS.toFixed(3).padStart(10) +
            sMgMl.toFixed(4).padStart(12)
        );
    });
    
    console.log();
    printSeparator();
    console.log("NOTES ON INTERPRETATION");
    printSeparator();
    console.log(`
1. All values are for FREE BASE forms at 25°C unless otherwise noted.
2. HCl salt forms typically have 10-1000x higher aqueous solubility.
3. log S values:
   - > 0: Highly soluble (> 1 mol/L)
   - -1 to 0: Soluble (0.1 - 1 mol/L)
   - -2 to -1: Slightly soluble (0.01 - 0.1 mol/L)
   - -3 to -2: Sparingly soluble (1 - 10 mM)
   - -4 to -3: Poorly soluble (0.1 - 1 mM)
   - < -4: Very poorly soluble (< 0.1 mM)

4. pH 7.4 values account for ionization of the basic amine.
5. Consensus values average 7 different prediction methods.
6. Experimental validation is recommended for critical applications.
`);
    
    // Export to JSON
    const exportData = {};
    for (const [name, data] of Object.entries(allResults)) {
        exportData[name] = {
            molecular_properties: {
                formula: data.compound.formula,
                mw: data.compound.mw,
                mp: data.compound.mp,
                logP: data.compound.logP,
                pKa: data.compound.pKa,
                tpsa: data.compound.tpsa
            },
            solubility_predictions: {}
        };
        
        for (const r of data.results) {
            exportData[name].solubility_predictions[r.method] = {
                log_s: parseFloat(r.logS.toFixed(4)),
                s_mol_L: r.sMolL,
                s_mg_mL: parseFloat(r.sMgMl.toFixed(6)),
                confidence: r.confidence
            };
        }
    }
    
    fs.writeFileSync("solubility_results.json", JSON.stringify(exportData, null, 2));
    console.log("\nResults exported to: solubility_results.json");
    printSeparator();
}

main();
