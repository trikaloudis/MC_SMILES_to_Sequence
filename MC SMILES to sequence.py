import streamlit as st
from rdkit import Chem

# --- The Logic Function (Same as before, just returns data) ---
def analyze_microcystin_advanced(smiles):
    mol = Chem.MolFromSmiles(smiles)
    
    if not mol:
        return None, "Error: Invalid SMILES string."

    patterns = {
        "Arg": "[NH]C(=[NH])[NH2]",
        "Leu": "[CH2]C(C)C",
        "Trp": "c1[nH]c2ccccc2c1",
        "Met": "[CH2]S[CH3]",
        "Hty_specific": "[CH2][CH2]c1ccc(O)cc1",
        "Hph_specific": "[CH2][CH2]c1ccccc1",
        "Tyr_core": "c1ccc(O)cc1",
        "Phe_core": "c1ccccc1",
    }

    counts = {}
    for name, smarts in patterns.items():
        substructure = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(substructure)
        counts[name] = len(matches)

    residues_found = []
    if counts["Arg"] > 0: residues_found.extend(["Arg"] * counts["Arg"])
    if counts["Leu"] > 0: residues_found.extend(["Leu"] * counts["Leu"])
    if counts["Trp"] > 0: residues_found.extend(["Trp"] * counts["Trp"])
    if counts["Met"] > 0: residues_found.extend(["Met"] * counts["Met"])

    hty_count = counts["Hty_specific"]
    if hty_count > 0: residues_found.extend(["Hty"] * hty_count)
    
    real_tyr_count = counts["Tyr_core"] - hty_count
    if real_tyr_count > 0: residues_found.extend(["Tyr"] * real_tyr_count)

    hph_count = counts["Hph_specific"]
    if hph_count > 0: residues_found.extend(["Hph"] * hph_count)

    # Subtract 1 for Adda phenyl ring
    real_phe_count = counts["Phe_core"] - 1 - hph_count
    if real_phe_count > 0: residues_found.extend(["Phe"] * real_phe_count)

    # Naming Logic
    X_residue = "?"
    Z_residue = "?"
    variant = "Unknown / Rare"

    if len(residues_found) >= 2:
        r1, r2 = residues_found[0], residues_found[1]
        code_map = {"Leu": "L", "Arg": "R", "Tyr": "Y", "Phe": "F", "Ala": "A", "Trp": "W", "Met": "M"}
        c1 = code_map.get(r1, "?")
        c2 = code_map.get(r2, "?")
        
        if "L" in [c1, c2] and "R" in [c1, c2]: variant = "MC-LR"
        elif "R" in [c1, c2] and c1==c2:         variant = "MC-RR"
        elif "Y" in [c1, c2] and "R" in [c1, c2]: variant = "MC-YR"
        elif "L" in [c1, c2] and "F" in [c1, c2]: variant = "MC-LF"
        elif "L" in [c1, c2] and "W" in [c1, c2]: variant = "MC-LW"
        else: variant = f"MC-{c1}{c2}"
        
        X_residue = r1
        Z_residue = r2
    elif "Leu" in residues_found and len(residues_found) == 1:
        variant = "MC-LA (Probable)"
        X_residue = "Leu"
        Z_residue = "Ala"

    return {
        "variant": variant,
        "sequence": f"Cyclo(D-Ala - {X_residue} - D-MeAsp - {Z_residue} - Adda - D-Glu - Mdha)",
        "residues": residues_found,
        "counts": counts
    }, None

# --- Streamlit UI Layout ---
st.set_page_config(page_title="Microcystin Analyzer", page_icon="🧬")

st.title("🧬 Microcystin SMILES Analyzer")
st.markdown("Paste a SMILES string below to identify the **Microcystin Variant** and **Amino Acid Sequence**.")

# Input area
user_smiles = st.text_area("Enter SMILES String:", height=100)

if st.button("Analyze Structure"):
    if user_smiles:
        result, error = analyze_microcystin_advanced(user_smiles)
        
        if error:
            st.error(error)
        else:
            # Display Results
            st.success(f"**Identified Variant:** {result['variant']}")
            
            st.subheader("Amino Acid Sequence")
            st.code(result['sequence'], language="text")
            
            st.subheader("Detailed Analysis")
            col1, col2 = st.columns(2)
            with col1:
                st.write("**Residues Detected:**")
                st.write(result['residues'])
            with col2:
                st.write("**Debug Counts:**")
                st.json(result['counts'])
                
            # Render Structure Image
            mol = Chem.MolFromSmiles(user_smiles)
            st.image(Chem.Draw.MolToImage(mol, size=(600, 400)), caption="Molecular Structure")
    else:
        st.warning("Please enter a SMILES string.")

st.sidebar.markdown("### Examples")
st.sidebar.code("MC-LR (Leu-Arg)", language="text")
st.sidebar.code("MC-RR (Arg-Arg)", language="text")