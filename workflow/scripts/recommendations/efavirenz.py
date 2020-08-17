def _recommendation(pt: str) -> str:
    recommendation = ""
    s1 = "cyp2b6"
    if s1 in pt:
        if pt[s1] == "ultrarapid_metabolizer":
            recommendation = "No recommendation"
        if pt[s1] == "rapid_metabolizer":
            recommendation = "No recommendation"
        if pt[s1] == "normal_metabolizer":
            recommendation = "No recommendation"
        elif pt[s1] == "intermediate_metabolizer":
            recommendation = "This patient is predicted to be a CYP2B6 intermediate metabolizer and may be at an increased risk of adverse events due to higher plasma concentrations of efavirenz. Consider initiating efavirenz with a decreased dose of 400 mg/day. Please consult a clinical pharmacistb for more information."
        elif pt[s1] == "poor_metabolizer":
            recommendation = "This patient is predicted to be a CYP2B6 poor metabolizer and may be at an increased risk of adverse events due to higher plasma concentrations of efavirenz. Consider initiating efavirenz with a decreased dose of 200 or 400 mg/day. Please consult a clinical pharmacistb for more information."
        else:
            recommendation = "No recommendation"
    elif s1 not in pt: 
        recommendation = "No recommendation"
    return recommendation