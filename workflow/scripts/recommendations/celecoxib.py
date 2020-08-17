def _recommendation(pt: str) -> str:
    recommendation = ""
    s1 = "cyp2c9"
    if s1 in pt:
        if pt[s1] == "ultrarapid_metabolizer":
            recommendation = "No recommendation"
        if pt[s1] == "rapid_metabolizer":
            recommendation = "No recommendation"
        if pt[s1] == "normal_metabolizer":
            recommendation = "No recommendation"
        elif pt[s1] == "intermediate_metabolizer":
            recommendation = "Based on the genotype result, this patient is predicted to be a CYP2C9 Intermediate Metabolizer with an activity score of 1 and is expected to have  higher plasma concentrations of celecoxib which may increase the risk of toxicities. Initiate therapy with lowest recommended starting dose. Titrate dose to clinical effect or maximum recommended dose with caution. In accordance with the prescribing information, use the lowest effective dosage for shortest duration consistent with individual patient treatment goals. Carefully monitor adverse events such as blood pressure and kidney function during course of therapy. Please consult a clinical pharmacist for more information."
        elif pt[s1] == "poor_metabolizer":
            recommendation = "Based on the genotype result, this patient is predicted to be a CYP2C9 Poor Metabolizer and is expected to have  higher plasma concentrations of celecoxib which may increase the risk of toxicities.  Initiate therapy with 25-50% of the lowest recommended starting dose. Titrate dose to clinical effect or 25-50% of the maximum recommended dose with caution. In accordance with the prescribing information, use the lowest effective dosage for shortest duration consistent with individual patient treatment goals. Dose titration should not occur until after steady state is reached (at least 8 days for celecoxib).  Carefully monitor adverse events such as blood pressure and kidney function during course of therapy. Alternatively, consider an alternate therapy not metabolized by CYP2C9 or not significantly impacted by CYP2C9 genetic variants in vivo. Please consult a clinical pharmacist for more information."
        else:
            recommendation = "No recommendation"
    elif s1 not in pt: 
        recommendation = "No recommendation"
    return recommendation