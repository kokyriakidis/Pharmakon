def Desipramine_recommendation(pt: str) -> str:
    recommendation = ""
    s1 = "cyp2d6"
    if s1 in pt:
        if pt[s1] == "ultrarapid_metabolizer":
            recommendation = "Based on the genotype result, this patient is predicted to be a CYP2D6 ultrarapid metabolizer and have lower than expected plasma concentrations resulting in an increased probability of pharmacotherapy failure. Consider alternative drug not metabolized by CYP2D6.  If Desipramine is warranted utilize therapeutic drug monitoring to guide dose adjustments.  Please consult a clinical pharmacist for more information."
        if pt[s1] == "normal_metabolizer":
            recommendation = "No recommendation"
        elif pt[s1] == "intermediate_metabolizer":
            recommendation = "Based on the genotype result, this patient is predicted to be a CYP2D6 intermediate metabolizer and higher than expected plasma concentrations of Desipramine resulting in an increased the probability of adverse reactions. Consider 25% reduction of recommended starting dose. Utilize therapeutic drug monitoring to guide dose adjustments. Please consult a clinical pharmacist for more information."
        elif pt[s1] == "poor_metabolizer":
            recommendation = "Based on the genotype result, this patient is predicted to be a CYP2D6 poor metabolizer and higher than expected plasma concentrations of Desipramine resulting in an increased the probability of adverse reactions. Consider 50% reduction of recommended starting dose. Utilize therapeutic drug monitoring to guide dose adjustments. Please consult a clinical pharmacist for more information."
        else:
            recommendation = "No recommendation"
    elif s1 not in pt: 
        recommendation = "No recommendation"
    return recommendation