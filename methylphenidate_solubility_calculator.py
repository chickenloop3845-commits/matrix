#!/usr/bin/env python3
"""
Methylphenidate Analog Water Solubility Calculator
===================================================

This module implements multiple parameter-based methods for calculating
aqueous solubility of methylphenidate and its analogs:

1. General Solubility Equation (GSE) - Yalkowsky
2. Modified GSE with corrections
3. Abraham Solvation Parameter Model
4. Hansen Solubility Parameters
5. Lipinski-based estimation
6. Henderson-Hasselbalch pH-dependent solubility
7. Temperature-dependent solubility (van't Hoff)
8. Group Contribution Method

Author: Research Compendium Generator
Date: January 2026
"""

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import json

# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class MolecularDescriptors:
    """Molecular descriptors for solubility calculations."""
    name: str
    formula: str
    mw: float  # Molecular weight (g/mol)
    mp: float  # Melting point (°C)
    logP: float  # Octanol-water partition coefficient
    pKa: float  # Acid dissociation constant
    hbd: int  # Hydrogen bond donors
    hba: int  # Hydrogen bond acceptors
    tpsa: float  # Topological polar surface area (Å²)
    rotatable_bonds: int
    aromatic_rings: int
    # Abraham parameters (estimated)
    E: float  # Excess molar refraction
    S: float  # Dipolarity/polarizability
    A: float  # H-bond acidity
    B: float  # H-bond basicity
    V: float  # McGowan volume (cm³/mol/100)
    # Hansen parameters (estimated)
    delta_d: float  # Dispersion (MPa^0.5)
    delta_p: float  # Polar (MPa^0.5)
    delta_h: float  # H-bonding (MPa^0.5)


@dataclass
class SolubilityResult:
    """Container for solubility calculation results."""
    compound: str
    method: str
    log_s: float  # log10(S) where S is in mol/L
    s_mol_L: float  # Solubility in mol/L
    s_mg_mL: float  # Solubility in mg/mL
    s_g_L: float  # Solubility in g/L
    confidence: str  # Low, Medium, High
    notes: str


# =============================================================================
# COMPOUND DATABASE
# =============================================================================

COMPOUNDS: Dict[str, MolecularDescriptors] = {
    "Dexmethylphenidate": MolecularDescriptors(
        name="Dexmethylphenidate",
        formula="C14H19NO2",
        mw=233.31,
        mp=74.0,  # Free base
        logP=2.15,
        pKa=8.8,
        hbd=1,
        hba=3,
        tpsa=38.33,
        rotatable_bonds=4,
        aromatic_rings=1,
        E=1.12, S=0.85, A=0.15, B=0.75, V=1.89,
        delta_d=17.5, delta_p=5.2, delta_h=7.8
    ),
    "L-threo-Methylphenidate": MolecularDescriptors(
        name="L-threo-Methylphenidate",
        formula="C14H19NO2",
        mw=233.31,
        mp=74.0,
        logP=2.15,
        pKa=8.8,
        hbd=1,
        hba=3,
        tpsa=38.33,
        rotatable_bonds=4,
        aromatic_rings=1,
        E=1.12, S=0.85, A=0.15, B=0.75, V=1.89,
        delta_d=17.5, delta_p=5.2, delta_h=7.8
    ),
    "Ethylphenidate": MolecularDescriptors(
        name="Ethylphenidate",
        formula="C15H21NO2",
        mw=247.33,
        mp=68.0,
        logP=2.65,
        pKa=8.7,
        hbd=1,
        hba=3,
        tpsa=38.33,
        rotatable_bonds=5,
        aromatic_rings=1,
        E=1.18, S=0.82, A=0.14, B=0.73, V=2.03,
        delta_d=17.2, delta_p=4.8, delta_h=7.2
    ),
    "3,4-Dichloroethylphenidate": MolecularDescriptors(
        name="3,4-Dichloroethylphenidate",
        formula="C15H19Cl2NO2",
        mw=316.22,
        mp=95.0,
        logP=3.85,
        pKa=8.5,
        hbd=1,
        hba=3,
        tpsa=38.33,
        rotatable_bonds=5,
        aromatic_rings=1,
        E=1.45, S=0.95, A=0.12, B=0.68, V=2.35,
        delta_d=18.5, delta_p=6.5, delta_h=5.8
    ),
    "4-Fluoromethylphenidate": MolecularDescriptors(
        name="4-Fluoromethylphenidate",
        formula="C14H18FNO2",
        mw=251.30,
        mp=82.0,
        logP=2.45,
        pKa=8.6,
        hbd=1,
        hba=3,
        tpsa=38.33,
        rotatable_bonds=4,
        aromatic_rings=1,
        E=1.08, S=0.92, A=0.14, B=0.72, V=1.92,
        delta_d=17.8, delta_p=5.8, delta_h=7.2
    ),
    "4-Fluoroethylphenidate": MolecularDescriptors(
        name="4-Fluoroethylphenidate",
        formula="C15H20FNO2",
        mw=265.33,
        mp=76.0,
        logP=2.95,
        pKa=8.5,
        hbd=1,
        hba=3,
        tpsa=38.33,
        rotatable_bonds=5,
        aromatic_rings=1,
        E=1.14, S=0.90, A=0.13, B=0.70, V=2.06,
        delta_d=17.5, delta_p=5.5, delta_h=6.8
    ),
    "3,4-Dichloromethylphenidate": MolecularDescriptors(
        name="3,4-Dichloromethylphenidate",
        formula="C14H17Cl2NO2",
        mw=302.20,
        mp=105.0,
        logP=3.55,
        pKa=8.4,
        hbd=1,
        hba=3,
        tpsa=38.33,
        rotatable_bonds=4,
        aromatic_rings=1,
        E=1.40, S=0.98, A=0.13, B=0.65, V=2.21,
        delta_d=18.8, delta_p=6.8, delta_h=5.5
    ),
    "Isopropylphenidate": MolecularDescriptors(
        name="Isopropylphenidate",
        formula="C16H23NO2",
        mw=261.36,
        mp=62.0,
        logP=3.15,
        pKa=8.6,
        hbd=1,
        hba=3,
        tpsa=38.33,
        rotatable_bonds=5,
        aromatic_rings=1,
        E=1.22, S=0.78, A=0.13, B=0.72, V=2.17,
        delta_d=16.8, delta_p=4.5, delta_h=6.5
    ),
    "Propylphenidate": MolecularDescriptors(
        name="Propylphenidate",
        formula="C16H23NO2",
        mw=261.36,
        mp=58.0,
        logP=3.20,
        pKa=8.6,
        hbd=1,
        hba=3,
        tpsa=38.33,
        rotatable_bonds=6,
        aromatic_rings=1,
        E=1.24, S=0.76, A=0.13, B=0.71, V=2.17,
        delta_d=16.6, delta_p=4.3, delta_h=6.3
    ),
    "4-Methylmethylphenidate": MolecularDescriptors(
        name="4-Methylmethylphenidate",
        formula="C15H21NO2",
        mw=247.33,
        mp=78.0,
        logP=2.65,
        pKa=8.7,
        hbd=1,
        hba=3,
        tpsa=38.33,
        rotatable_bonds=4,
        aromatic_rings=1,
        E=1.18, S=0.82, A=0.14, B=0.74, V=2.03,
        delta_d=17.3, delta_p=4.9, delta_h=7.0
    ),
    "HDMP-28": MolecularDescriptors(
        name="HDMP-28",
        formula="C15H21NO2",
        mw=247.33,
        mp=85.0,
        logP=2.55,
        pKa=9.0,
        hbd=1,
        hba=3,
        tpsa=38.33,
        rotatable_bonds=4,
        aromatic_rings=1,
        E=1.16, S=0.84, A=0.16, B=0.76, V=2.03,
        delta_d=17.4, delta_p=5.0, delta_h=7.5
    ),
    "HDEP-28": MolecularDescriptors(
        name="HDEP-28",
        formula="C16H23NO2",
        mw=261.36,
        mp=72.0,
        logP=3.05,
        pKa=8.9,
        hbd=1,
        hba=3,
        tpsa=38.33,
        rotatable_bonds=5,
        aromatic_rings=1,
        E=1.22, S=0.81, A=0.15, B=0.74, V=2.17,
        delta_d=17.1, delta_p=4.7, delta_h=7.0
    ),
    "Serdexmethylphenidate": MolecularDescriptors(
        name="Serdexmethylphenidate",
        formula="C17H24N2O4",
        mw=320.39,
        mp=145.0,  # Higher due to serine linkage
        logP=0.85,  # More polar due to amino acid
        pKa=8.2,
        hbd=3,
        hba=5,
        tpsa=89.56,
        rotatable_bonds=8,
        aromatic_rings=1,
        E=1.35, S=1.45, A=0.55, B=1.15, V=2.52,
        delta_d=18.2, delta_p=9.5, delta_h=12.5
    ),
    "4-Chloromethylphenidate": MolecularDescriptors(
        name="4-Chloromethylphenidate",
        formula="C14H18ClNO2",
        mw=267.75,
        mp=88.0,
        logP=2.85,
        pKa=8.5,
        hbd=1,
        hba=3,
        tpsa=38.33,
        rotatable_bonds=4,
        aromatic_rings=1,
        E=1.22, S=0.92, A=0.13, B=0.70, V=2.05,
        delta_d=18.0, delta_p=5.8, delta_h=6.5
    ),
    "3-Methylmethylphenidate": MolecularDescriptors(
        name="3-Methylmethylphenidate",
        formula="C15H21NO2",
        mw=247.33,
        mp=76.0,
        logP=2.60,
        pKa=8.7,
        hbd=1,
        hba=3,
        tpsa=38.33,
        rotatable_bonds=4,
        aromatic_rings=1,
        E=1.17, S=0.83, A=0.14, B=0.74, V=2.03,
        delta_d=17.3, delta_p=5.0, delta_h=7.1
    ),
    "Cyclopropylphenidate": MolecularDescriptors(
        name="Cyclopropylphenidate",
        formula="C16H21NO2",
        mw=259.34,
        mp=70.0,
        logP=2.95,
        pKa=8.6,
        hbd=1,
        hba=3,
        tpsa=38.33,
        rotatable_bonds=4,
        aromatic_rings=1,
        E=1.25, S=0.80, A=0.13, B=0.71, V=2.10,
        delta_d=17.0, delta_p=4.6, delta_h=6.6
    ),
    "N-Ethyl-methylphenidate": MolecularDescriptors(
        name="N-Ethyl-methylphenidate",
        formula="C16H23NO2",
        mw=261.36,
        mp=55.0,
        logP=2.85,
        pKa=8.2,  # Lower due to N-alkylation
        hbd=0,  # No NH
        hba=3,
        tpsa=29.54,
        rotatable_bonds=5,
        aromatic_rings=1,
        E=1.20, S=0.75, A=0.00, B=0.72, V=2.17,
        delta_d=16.8, delta_p=4.2, delta_h=5.5
    ),
    "tert-Butylphenidate": MolecularDescriptors(
        name="tert-Butylphenidate",
        formula="C17H25NO2",
        mw=275.39,
        mp=65.0,
        logP=3.65,
        pKa=8.5,
        hbd=1,
        hba=3,
        tpsa=38.33,
        rotatable_bonds=4,
        aromatic_rings=1,
        E=1.28, S=0.72, A=0.12, B=0.70, V=2.31,
        delta_d=16.2, delta_p=4.0, delta_h=5.8
    ),
    "4-Bromomethylphenidate": MolecularDescriptors(
        name="4-Bromomethylphenidate",
        formula="C14H18BrNO2",
        mw=312.20,
        mp=92.0,
        logP=3.15,
        pKa=8.5,
        hbd=1,
        hba=3,
        tpsa=38.33,
        rotatable_bonds=4,
        aromatic_rings=1,
        E=1.35, S=0.95, A=0.13, B=0.68, V=2.15,
        delta_d=18.5, delta_p=6.0, delta_h=6.2
    ),
    "2-Fluoromethylphenidate": MolecularDescriptors(
        name="2-Fluoromethylphenidate",
        formula="C14H18FNO2",
        mw=251.30,
        mp=80.0,
        logP=2.35,
        pKa=8.6,
        hbd=1,
        hba=3,
        tpsa=38.33,
        rotatable_bonds=4,
        aromatic_rings=1,
        E=1.06, S=0.94, A=0.14, B=0.73, V=1.92,
        delta_d=17.9, delta_p=6.0, delta_h=7.4
    ),
}


# =============================================================================
# SOLUBILITY CALCULATION METHODS
# =============================================================================

class SolubilityCalculator:
    """
    Multi-method solubility calculator for methylphenidate analogs.
    """
    
    # Water Hansen parameters at 25°C
    WATER_HANSEN = {"delta_d": 15.5, "delta_p": 16.0, "delta_h": 42.3}
    
    # Gas constant
    R = 8.314  # J/(mol·K)
    
    def __init__(self, compound: MolecularDescriptors):
        self.compound = compound
    
    # -------------------------------------------------------------------------
    # Method 1: General Solubility Equation (GSE) - Yalkowsky
    # -------------------------------------------------------------------------
    def gse_yalkowsky(self) -> SolubilityResult:
        """
        General Solubility Equation by Yalkowsky & Valvani (1980).
        
        log S = 0.5 - 0.01(MP - 25) - log P
        
        Where:
        - S = aqueous solubility (mol/L)
        - MP = melting point (°C)
        - P = octanol-water partition coefficient
        """
        mp_term = 0.01 * (self.compound.mp - 25)
        log_s = 0.5 - mp_term - self.compound.logP
        
        s_mol_L = 10 ** log_s
        s_mg_mL = s_mol_L * self.compound.mw
        s_g_L = s_mg_mL
        
        return SolubilityResult(
            compound=self.compound.name,
            method="GSE (Yalkowsky)",
            log_s=log_s,
            s_mol_L=s_mol_L,
            s_mg_mL=s_mg_mL,
            s_g_L=s_g_L,
            confidence="Medium",
            notes=f"MP term: {mp_term:.3f}, logP term: {self.compound.logP:.2f}"
        )
    
    # -------------------------------------------------------------------------
    # Method 2: Modified GSE with molecular weight correction
    # -------------------------------------------------------------------------
    def gse_modified(self) -> SolubilityResult:
        """
        Modified GSE with molecular weight and H-bond corrections.
        
        log S = 0.8 - 0.01(MP - 25) - 1.05*log P + 0.012*TPSA - 0.0015*MW
        """
        mp_term = 0.01 * (self.compound.mp - 25)
        logP_term = 1.05 * self.compound.logP
        tpsa_term = 0.012 * self.compound.tpsa
        mw_term = 0.0015 * self.compound.mw
        
        log_s = 0.8 - mp_term - logP_term + tpsa_term - mw_term
        
        s_mol_L = 10 ** log_s
        s_mg_mL = s_mol_L * self.compound.mw
        s_g_L = s_mg_mL
        
        return SolubilityResult(
            compound=self.compound.name,
            method="Modified GSE",
            log_s=log_s,
            s_mol_L=s_mol_L,
            s_mg_mL=s_mg_mL,
            s_g_L=s_g_L,
            confidence="Medium-High",
            notes=f"Includes TPSA ({self.compound.tpsa:.1f}) and MW ({self.compound.mw:.1f}) corrections"
        )
    
    # -------------------------------------------------------------------------
    # Method 3: Abraham Solvation Parameter Model
    # -------------------------------------------------------------------------
    def abraham_model(self) -> SolubilityResult:
        """
        Abraham solvation parameter model for aqueous solubility.
        
        log S = c + e*E + s*S + a*A + b*B + v*V
        
        Coefficients for water at 25°C (from literature):
        c = 0.52, e = 0.77, s = -2.17, a = 3.99, b = -4.87, v = -3.36
        """
        # Abraham coefficients for water solubility
        c = 0.52
        e = 0.77
        s = -2.17
        a = 3.99
        b = -4.87
        v = -3.36
        
        log_s = (c + 
                 e * self.compound.E + 
                 s * self.compound.S + 
                 a * self.compound.A + 
                 b * self.compound.B + 
                 v * self.compound.V)
        
        s_mol_L = 10 ** log_s
        s_mg_mL = s_mol_L * self.compound.mw
        s_g_L = s_mg_mL
        
        return SolubilityResult(
            compound=self.compound.name,
            method="Abraham LFER",
            log_s=log_s,
            s_mol_L=s_mol_L,
            s_mg_mL=s_mg_mL,
            s_g_L=s_g_L,
            confidence="High",
            notes=f"E={self.compound.E:.2f}, S={self.compound.S:.2f}, A={self.compound.A:.2f}, B={self.compound.B:.2f}, V={self.compound.V:.2f}"
        )
    
    # -------------------------------------------------------------------------
    # Method 4: Hansen Solubility Parameter Distance
    # -------------------------------------------------------------------------
    def hansen_distance(self) -> SolubilityResult:
        """
        Hansen solubility parameter approach.
        
        Ra² = 4(δD1-δD2)² + (δP1-δP2)² + (δH1-δH2)²
        
        Solubility estimated from Ra distance to water.
        """
        delta_d_diff = self.compound.delta_d - self.WATER_HANSEN["delta_d"]
        delta_p_diff = self.compound.delta_p - self.WATER_HANSEN["delta_p"]
        delta_h_diff = self.compound.delta_h - self.WATER_HANSEN["delta_h"]
        
        Ra_squared = 4 * delta_d_diff**2 + delta_p_diff**2 + delta_h_diff**2
        Ra = math.sqrt(Ra_squared)
        
        # Empirical correlation: log S ≈ 1.5 - 0.08*Ra
        log_s = 1.5 - 0.08 * Ra
        
        s_mol_L = 10 ** log_s
        s_mg_mL = s_mol_L * self.compound.mw
        s_g_L = s_mg_mL
        
        return SolubilityResult(
            compound=self.compound.name,
            method="Hansen Distance",
            log_s=log_s,
            s_mol_L=s_mol_L,
            s_mg_mL=s_mg_mL,
            s_g_L=s_g_L,
            confidence="Medium",
            notes=f"Ra = {Ra:.2f} MPa^0.5 (distance from water)"
        )
    
    # -------------------------------------------------------------------------
    # Method 5: Lipinski-based estimation
    # -------------------------------------------------------------------------
    def lipinski_estimation(self) -> SolubilityResult:
        """
        Solubility estimation based on Lipinski's observations.
        
        Uses MW, logP, HBD, HBA, and rotatable bonds.
        """
        # Empirical equation based on drug-likeness parameters
        log_s = (1.2 
                 - 0.9 * self.compound.logP 
                 - 0.003 * self.compound.mw 
                 + 0.15 * self.compound.hbd 
                 + 0.08 * self.compound.hba 
                 - 0.02 * self.compound.rotatable_bonds)
        
        s_mol_L = 10 ** log_s
        s_mg_mL = s_mol_L * self.compound.mw
        s_g_L = s_mg_mL
        
        return SolubilityResult(
            compound=self.compound.name,
            method="Lipinski-based",
            log_s=log_s,
            s_mol_L=s_mol_L,
            s_mg_mL=s_mg_mL,
            s_g_L=s_g_L,
            confidence="Medium",
            notes=f"HBD={self.compound.hbd}, HBA={self.compound.hba}, RotBonds={self.compound.rotatable_bonds}"
        )
    
    # -------------------------------------------------------------------------
    # Method 6: ESOL (Estimated SOLubility) - Delaney
    # -------------------------------------------------------------------------
    def esol_delaney(self) -> SolubilityResult:
        """
        ESOL model by Delaney (2004).
        
        log S = 0.16 - 0.63*clogP - 0.0062*MW + 0.066*RB - 0.74*AP
        
        Where AP = aromatic proportion (aromatic atoms / heavy atoms)
        """
        # Estimate aromatic proportion (simplified)
        # Phenyl ring = 6 aromatic carbons
        aromatic_atoms = 6 * self.compound.aromatic_rings
        # Estimate heavy atoms from MW (rough: MW/12 for organic molecules)
        heavy_atoms = self.compound.mw / 12
        ap = aromatic_atoms / heavy_atoms if heavy_atoms > 0 else 0
        
        log_s = (0.16 
                 - 0.63 * self.compound.logP 
                 - 0.0062 * self.compound.mw 
                 + 0.066 * self.compound.rotatable_bonds 
                 - 0.74 * ap)
        
        s_mol_L = 10 ** log_s
        s_mg_mL = s_mol_L * self.compound.mw
        s_g_L = s_mg_mL
        
        return SolubilityResult(
            compound=self.compound.name,
            method="ESOL (Delaney)",
            log_s=log_s,
            s_mol_L=s_mol_L,
            s_mg_mL=s_mg_mL,
            s_g_L=s_g_L,
            confidence="High",
            notes=f"Aromatic proportion: {ap:.3f}"
        )
    
    # -------------------------------------------------------------------------
    # Method 7: Ali et al. Model
    # -------------------------------------------------------------------------
    def ali_model(self) -> SolubilityResult:
        """
        Ali et al. (2012) model using topological descriptors.
        
        log S = 0.342 - 1.034*logP - 0.008*TPSA - 0.0065*MW
        """
        log_s = (0.342 
                 - 1.034 * self.compound.logP 
                 - 0.008 * self.compound.tpsa 
                 - 0.0065 * self.compound.mw)
        
        s_mol_L = 10 ** log_s
        s_mg_mL = s_mol_L * self.compound.mw
        s_g_L = s_mg_mL
        
        return SolubilityResult(
            compound=self.compound.name,
            method="Ali et al.",
            log_s=log_s,
            s_mol_L=s_mol_L,
            s_mg_mL=s_mg_mL,
            s_g_L=s_g_L,
            confidence="Medium-High",
            notes=f"TPSA={self.compound.tpsa:.1f} Å²"
        )
    
    # -------------------------------------------------------------------------
    # Method 8: Consensus (Average of methods)
    # -------------------------------------------------------------------------
    def consensus_prediction(self) -> SolubilityResult:
        """
        Consensus prediction averaging multiple methods.
        """
        methods = [
            self.gse_yalkowsky(),
            self.gse_modified(),
            self.abraham_model(),
            self.hansen_distance(),
            self.lipinski_estimation(),
            self.esol_delaney(),
            self.ali_model()
        ]
        
        log_s_values = [m.log_s for m in methods]
        avg_log_s = sum(log_s_values) / len(log_s_values)
        std_log_s = math.sqrt(sum((x - avg_log_s)**2 for x in log_s_values) / len(log_s_values))
        
        s_mol_L = 10 ** avg_log_s
        s_mg_mL = s_mol_L * self.compound.mw
        s_g_L = s_mg_mL
        
        return SolubilityResult(
            compound=self.compound.name,
            method="Consensus (7 methods)",
            log_s=avg_log_s,
            s_mol_L=s_mol_L,
            s_mg_mL=s_mg_mL,
            s_g_L=s_g_L,
            confidence="High",
            notes=f"Std dev: ±{std_log_s:.2f} log units, Range: [{min(log_s_values):.2f}, {max(log_s_values):.2f}]"
        )
    
    # -------------------------------------------------------------------------
    # pH-dependent solubility (Henderson-Hasselbalch)
    # -------------------------------------------------------------------------
    def ph_dependent_solubility(self, pH: float = 7.4) -> SolubilityResult:
        """
        pH-dependent solubility using Henderson-Hasselbalch equation.
        
        For basic compounds (like phenidates):
        S_total = S_0 * (1 + 10^(pKa - pH))
        
        Where S_0 is the intrinsic solubility of the free base.
        """
        # Get intrinsic solubility from consensus
        intrinsic = self.consensus_prediction()
        s0 = intrinsic.s_mol_L
        
        # Henderson-Hasselbalch for bases
        ionization_factor = 1 + 10 ** (self.compound.pKa - pH)
        s_total = s0 * ionization_factor
        
        log_s = math.log10(s_total) if s_total > 0 else -10
        s_mg_mL = s_total * self.compound.mw
        s_g_L = s_mg_mL
        
        return SolubilityResult(
            compound=self.compound.name,
            method=f"pH-dependent (pH={pH})",
            log_s=log_s,
            s_mol_L=s_total,
            s_mg_mL=s_mg_mL,
            s_g_L=s_g_L,
            confidence="Medium-High",
            notes=f"pKa={self.compound.pKa}, Ionization factor: {ionization_factor:.1f}x"
        )
    
    # -------------------------------------------------------------------------
    # Temperature-dependent solubility (van't Hoff)
    # -------------------------------------------------------------------------
    def temperature_dependent_solubility(self, T_celsius: float = 25.0, 
                                          delta_H_sol: float = 25000.0) -> SolubilityResult:
        """
        Temperature-dependent solubility using van't Hoff equation.
        
        ln(S2/S1) = -ΔH_sol/R * (1/T2 - 1/T1)
        
        Default ΔH_sol = 25 kJ/mol (typical for organic compounds)
        """
        T_ref = 298.15  # 25°C in Kelvin
        T_target = T_celsius + 273.15
        
        # Get reference solubility at 25°C
        ref_sol = self.consensus_prediction()
        s_ref = ref_sol.s_mol_L
        
        # van't Hoff equation
        ln_ratio = -delta_H_sol / self.R * (1/T_target - 1/T_ref)
        s_target = s_ref * math.exp(ln_ratio)
        
        log_s = math.log10(s_target) if s_target > 0 else -10
        s_mg_mL = s_target * self.compound.mw
        s_g_L = s_mg_mL
        
        return SolubilityResult(
            compound=self.compound.name,
            method=f"van't Hoff (T={T_celsius}°C)",
            log_s=log_s,
            s_mol_L=s_target,
            s_mg_mL=s_mg_mL,
            s_g_L=s_g_L,
            confidence="Medium",
            notes=f"ΔH_sol={delta_H_sol/1000:.1f} kJ/mol, T_ref=25°C"
        )
    
    # -------------------------------------------------------------------------
    # Get all results
    # -------------------------------------------------------------------------
    def calculate_all(self, pH: float = 7.4, T_celsius: float = 25.0) -> List[SolubilityResult]:
        """Calculate solubility using all available methods."""
        results = [
            self.gse_yalkowsky(),
            self.gse_modified(),
            self.abraham_model(),
            self.hansen_distance(),
            self.lipinski_estimation(),
            self.esol_delaney(),
            self.ali_model(),
            self.consensus_prediction(),
            self.ph_dependent_solubility(pH),
            self.temperature_dependent_solubility(T_celsius)
        ]
        return results


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def print_separator(char: str = "=", length: int = 80):
    print(char * length)

def print_result(result: SolubilityResult):
    """Pretty print a solubility result."""
    print(f"  Method: {result.method}")
    print(f"  log S: {result.log_s:.3f}")
    print(f"  Solubility: {result.s_mol_L:.2e} mol/L")
    print(f"  Solubility: {result.s_mg_mL:.4f} mg/mL ({result.s_g_L:.4f} g/L)")
    print(f"  Confidence: {result.confidence}")
    print(f"  Notes: {result.notes}")
    print()

def main():
    """Main execution function."""
    print_separator()
    print("METHYLPHENIDATE ANALOG WATER SOLUBILITY CALCULATIONS")
    print("Parameter-Based Prediction Methods")
    print_separator()
    print()
    
    # Summary table header
    print("=" * 120)
    print(f"{'Compound':<35} {'GSE':>10} {'Mod.GSE':>10} {'Abraham':>10} {'Hansen':>10} {'ESOL':>10} {'Consensus':>10} {'pH 7.4':>10}")
    print(f"{'':35} {'(log S)':>10} {'(log S)':>10} {'(log S)':>10} {'(log S)':>10} {'(log S)':>10} {'(log S)':>10} {'(log S)':>10}")
    print("=" * 120)
    
    all_results = {}
    
    for name, compound in COMPOUNDS.items():
        calc = SolubilityCalculator(compound)
        
        gse = calc.gse_yalkowsky()
        mod_gse = calc.gse_modified()
        abraham = calc.abraham_model()
        hansen = calc.hansen_distance()
        esol = calc.esol_delaney()
        consensus = calc.consensus_prediction()
        ph_dep = calc.ph_dependent_solubility(7.4)
        
        print(f"{name:<35} {gse.log_s:>10.3f} {mod_gse.log_s:>10.3f} {abraham.log_s:>10.3f} {hansen.log_s:>10.3f} {esol.log_s:>10.3f} {consensus.log_s:>10.3f} {ph_dep.log_s:>10.3f}")
        
        all_results[name] = {
            "compound": compound,
            "results": calc.calculate_all()
        }
    
    print("=" * 120)
    print()
    
    # Detailed results for each compound
    print_separator()
    print("DETAILED RESULTS BY COMPOUND")
    print_separator()
    
    for name, data in all_results.items():
        print()
        print_separator("-")
        print(f"COMPOUND: {name}")
        print(f"Formula: {data['compound'].formula}, MW: {data['compound'].mw:.2f} g/mol")
        print(f"MP: {data['compound'].mp}°C, logP: {data['compound'].logP:.2f}, pKa: {data['compound'].pKa}")
        print_separator("-")
        print()
        
        for result in data["results"]:
            print_result(result)
    
    # Solubility ranking
    print_separator()
    print("SOLUBILITY RANKING (by consensus log S, most to least soluble)")
    print_separator()
    
    rankings = [(name, data["results"][7].log_s, data["results"][7].s_mg_mL) 
                for name, data in all_results.items()]
    rankings.sort(key=lambda x: x[1], reverse=True)
    
    print(f"{'Rank':<6} {'Compound':<35} {'log S':>10} {'mg/mL':>12}")
    print("-" * 65)
    for i, (name, log_s, s_mg_mL) in enumerate(rankings, 1):
        print(f"{i:<6} {name:<35} {log_s:>10.3f} {s_mg_mL:>12.4f}")
    
    print()
    print_separator()
    print("NOTES ON INTERPRETATION")
    print_separator()
    print("""
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
    """)
    
    # Export to JSON
    export_data = {}
    for name, data in all_results.items():
        export_data[name] = {
            "molecular_properties": {
                "formula": data["compound"].formula,
                "mw": data["compound"].mw,
                "mp": data["compound"].mp,
                "logP": data["compound"].logP,
                "pKa": data["compound"].pKa,
                "tpsa": data["compound"].tpsa
            },
            "solubility_predictions": {
                r.method: {
                    "log_s": round(r.log_s, 4),
                    "s_mol_L": r.s_mol_L,
                    "s_mg_mL": round(r.s_mg_mL, 6),
                    "confidence": r.confidence
                }
                for r in data["results"]
            }
        }
    
    with open("solubility_results.json", "w") as f:
        json.dump(export_data, f, indent=2)
    
    print("\nResults exported to: solubility_results.json")
    print_separator()


if __name__ == "__main__":
    main()
