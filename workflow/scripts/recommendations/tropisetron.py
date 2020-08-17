def _recommendation(pt: str) -> str:
    recommendation = ""
    s1 = "cyp2d6"
    if s1 in pt:
        if pt[s1] == "ultrarapid_metabolizer":
            recommendation = "This patient is predicted to be a CYP2D6 ultrarapid metabolizer and may be at an increased risk of a poor response due to low plasma concentrations of tropisetron. Consider selecting an alternative medication not extensively metabolized by CYP2D6 (i.e., granisetron). Please consult a clinical pharmacist for more information."
        if pt[s1] == "normal_metabolizer":
            recommendation = "No recommendation"
        elif pt[s1] == "intermediate_metabolizer":
            recommendation = "No recommendation"
        elif pt[s1] == "poor_metabolizer":
            recommendation = "No recommendation"
        else:
            recommendation = "No recommendation"
    elif s1 not in pt: 
        recommendation = "No recommendation"
    return recommendation