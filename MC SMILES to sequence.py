import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
import re

# ==========================================
# CORE LOGIC: SMILES ANALYSIS
# ==========================================
def analyze_smiles_chemistry(smiles):
    """
    Analyzes SMILES to find X and Z residues using chemical substructure matching.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None, "Invalid SMILES"

    # SMARTS patterns for amino acid side chains
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

    # DEDUCTION LOGIC
    residues_found = []
    if counts["Arg"] > 0: residues_found.extend(["Arg"] * counts["Arg"])
    if counts["Leu"] > 0: residues_found.extend(["Leu"] * counts["Leu"])
    if counts["Trp"] > 0: residues_found.extend(["Trp"] * counts["Trp"])
    if counts["Met"] > 0: residues_found.extend(["Met"] * counts["Met"])

    # Homo-amino acid handling
    hty_count = counts["Hty_specific"]
    if hty_count > 0: residues_found.extend(["Hty"] * hty_count)
    
    real_tyr_count = counts["Tyr_core"] - hty_count
    if real_tyr_count > 0: residues_found.extend(["Tyr"] * real_tyr_count)

    hph_count = counts["Hph_specific"]
    if hph_count > 0: residues_found.extend(["Hph"] * hph_count)

    # Subtract 1 Phenyl ring for the Adda moiety (always present)
    real_phe_count = counts["Phe_core"] - 1 - hph_count
    if real_phe_count > 0: residues_found.extend(["Phe"] * real_phe_count)

    # Determine X and Z candidates
    X_cand, Z_cand = "?", "?"
    
    # Sort to standard order (Hydrophobic first usually)
    # Heuristic: If Leu is present, it's usually X. Arg is usually Z in LR.
    # If two Args, both are X/Z.
    residues_found.sort(key=lambda x: 0 if x == "Leu" else 1)
    
    if len(residues_found) >= 2:
        X_cand = residues_found[0]
        Z_cand = residues_found[1]
    elif len(residues_found) == 1 and "Leu" in residues_found:
        # Implicit LA case
        X_cand = "Leu"
        Z_cand = "Ala (inferred)"

    return {
        "residues": residues_found,
        "X_chem": X_cand,
        "Z_chem": Z_cand
    }, None

# ==========================================
# CORE LOGIC: NAME PARSING
# ==========================================
def parse_compound_name(name):
    """
    Extracts text-based modifications (e.g. [D-Asp3]) from the Name Column.
    """
    # 1. Extract Modification Tags (content inside brackets)
    mods = re.findall(r"\[(.*?)\]", name)
    
    # 2. Extract Variant Suffix (e.g. "-LR", "-RR", "-LF")
    # Looks for pattern "Microcystin-XY" or "MC-XY"
    variant_match = re.search(r"Microcystin-([A-Z]{2})|MC-([A-Z]{2})", name, re.IGNORECASE)
    
    variant_suffix = None
    if variant_match:
        # Get the group that matched
        variant_suffix = variant_match.group(1) or variant_match.group(2)
    
    return {
        "mods": mods,
        "variant_suffix": variant_suffix
    }

# ==========================================
# MASTER ROW PROCESSOR
# ==========================================
def process_microcystin_row(row):
    """
    Integrates Name info and SMILES info to build the final sequence and confidence.
    """
    # 1. Get Inputs
    # Try to find columns case-insensitively or fall back to indices
    col_map = {c.lower(): c for c in row.index}
    
    name = row.get(col_map.get("compoundname"), "")
    smiles = row.get(col_map.get("smiles"), "")
    
    # If SMILES is empty, check column J (index 9) if it exists
    if not smiles and len(row) > 9:
         smiles = row.iloc[9]

    if not isinstance(smiles, str) or not smiles.strip():
        return pd.Series(["Error", "Missing SMILES", "Low", "No SMILES provided"])

    # 2. Run Analyses
    chem_data, error = analyze_smiles_chemistry(smiles)
    if error:
         return pd.Series(["Error", "Invalid SMILES", "Low", "RDKit could not parse SMILES"])
         
    text_data = parse_compound_name(str(name))
    
    # 3. Build Sequence Template
    # Standard: Cyclo(D-Ala1 - X2 - D-MeAsp3 - Z4 - Adda5 - D-Glu6 - Mdha7)
    pos1 = "D-Ala"
    pos3 = "D-MeAsp"
    pos6 = "D-Glu"
    pos7 = "Mdha"
    
    X_final = chem_data["X_chem"]
    Z_final = chem_data["Z_chem"]
    
    # 4. Apply Name Modifications
    comments = []
    
    # Check for [D-Asp3] (Desmethyl)
    for mod in text_data["mods"]:
        if "Asp3" in mod:
            pos3 = "D-Asp"
            comments.append("Modification: Desmethyl-Asp at Pos 3")
        if "Dha7" in mod:
            pos7 = "Dha"
            comments.append("Modification: Dehydroalanine at Pos 7")
        if "Mser7" in mod:
            pos7 = "Mser"
            comments.append("Modification: Methylserine at Pos 7")
        if "Ser7" in mod:
            pos7 = "Ser"
            comments.append("Modification: Serine at Pos 7")
            
    # 5. Consistency Check (Confidence)
    confidence = "High"
    
    # Check if Name Variant matches SMILES Chemistry
    # Map letters to full names
    aa_map = {"L": "Leu", "R": "Arg", "Y": "Tyr", "F": "Phe", "W": "Trp", "M": "Met", "A": "Ala"}
    
    if text_data["variant_suffix"]:
        name_X = aa_map.get(text_data["variant_suffix"][0], "?")
        name_Z = aa_map.get(text_data["variant_suffix"][1], "?")
        
        # We check if our Chem findings match the Name
        # Note: Order in chem findings might be swapped, so we check set membership
        chem_set = {X_final, Z_final}
        name_set = {name_X, name_Z}
        
        # Simple intersection check
        # We are lenient: if at least one matches, or if order is just swapped
        if chem_set == name_set:
            comments.append("Perfect Match: Name and SMILES agree on X/Z.")
        elif len(chem_set.intersection(name_set)) > 0:
            confidence = "Medium"
            comments.append(f"Partial Match: Name suggests {name_X}/{name_Z}, SMILES found {X_final}/{Z_final}.")
        else:
            confidence = "Low"
            comments.append(f"Mismatch: Name implies {name_X}/{name_Z} but SMILES contains {X_final}/{Z_final}.")
    else:
        # No variant suffix found in name
        confidence = "Medium"
        comments.append("Variant not specified in name (e.g. -LR), relying on SMILES only.")

    sequence = f"Cyclo({pos1} - {X_final} - {pos3} - {Z_final} - Adda - {pos6} - {pos7})"
    
    return pd.Series([sequence, f"MC-{X_final}{Z_final} (SMILES)", confidence, "; ".join(comments)])


# ==========================================
# STREAMLIT UI
# ==========================================
st.set_page_config(page_title="CyanoMetDB Analyzer", page_icon="🧪", layout="wide")

st.title("🧪 CyanoMetDB Batch Processor")
st.markdown("""
**Upload your CSV file containing Microcystin data.** The tool will analyze the **CompoundName** (Col B) and **SMILES** (Col J) to generate the amino acid sequence.
""")

uploaded_file = st.file_uploader("Upload CSV File", type=["csv"])

if uploaded_file is not None:
    try:
        # Load Data
        df = pd.read_csv(uploaded_file)
        
        # Validate Columns (Loose check)
        # We need to ensure we can process it even if headers differ slightly
        st.write("Preview of Uploaded Data:")
        st.dataframe(df.head(3))
        
        if st.button("Run Analysis"):
            with st.spinner("Analyzing SMILES patterns..."):
                # Apply the processor
                # We expect new columns: Sequence, Detected_Variant, Confidence, Comments
                new_cols = df.apply(process_microcystin_row, axis=1)
                new_cols.columns = ["Proposed_Sequence", "SMILES_Variant", "Confidence", "Analysis_Notes"]
                
                # Combine original df with new results
                final_df = pd.concat([df, new_cols], axis=1)
                
                st.success("Analysis Complete!")
                
                # Show results with highlighting
                st.write("Results Preview:")
                st.dataframe(final_df.head(10))
                
                # CSV Download
                csv = final_df.to_csv(index=False).encode('utf-8')
                st.download_button(
                    label="📥 Download Processed Table",
                    data=csv,
                    file_name="microcystin_analysis_results.csv",
                    mime="text/csv",
                )
                
    except Exception as e:
        st.error(f"Error processing file: {e}")
