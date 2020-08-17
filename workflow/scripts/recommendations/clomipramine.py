def Clomipramine_recommendation(pt: str) -> str:
    recommendation = ""
    s1 = "cyp2d6"
    s2 = "cyp2c19"
    if s1 in pt and s2 in pt:
        if pt[s1] == "ultrarapid_metabolizer" and pt[s2] == "ultrarapid_metabolizer":
            recommendation = "This patient is predicated to be a CYP2D6 ultrarapid metabolizer and a CYP2C19 ultrarapid metabolizer. This patient may be at an increased risk of a sub-optimal response. Consider selecting an alternative drug.  If Clomipramine is warranted utilize therapeutic drug monitoring to guide dose adjustments. Please consult a clinical pharmacist for more information."
        elif pt[s1] == "ultrarapid_metabolizer" and pt[s2] == "rapid_metabolizer":
            recommendation = "This patient is predicated to be a CYP2D6 ultrarapid metabolizer and a CYP2C19 rapid metabolizer. This patient may be at an increased risk of a sub-optimal response. Consider selecting an alternative drug.  If Clomipramine is warranted utilize therapeutic drug monitoring to guide dose adjustments. Please consult a clinical pharmacist for more information."
        elif pt[s1] == "ultrarapid_metabolizer" and pt[s2] == "normal_metabolizer":
            recommendation = "This patient is predicted to be a CYP2D6 ultrarapid metabolizer and CYP2C19 normal metabolizer. Based on the CYP2D6 result, the patient may be at an increased risk of a sub-optimal response. Consider selecting an alternative drug not metabolized by CYP2D6. If Clomipramine is warranted utilize therapeutic drug monitoring to guide dose adjustments. Please consult a clinical pharmacist for more information."
        elif pt[s1] == "ultrarapid_metabolizer" and pt[s2] == "intermediate_metabolizer":
            recommendation = "This patient is predicated to be a CYP2D6 ultrarapid metabolizer and a CYP2C19 intermediate metabolizer. This patient may be at an increased risk of a sub-optimal response. Consider selecting an alternative drug.  If Clomipramine is warranted utilize therapeutic drug monitoring to guide dose adjustments. Please consult a clinical pharmacist for more information."
        elif pt[s1] == "ultrarapid_metabolizer" and pt[s2] == "poor_metabolizer":
            recommendation = "This patient is predicated to be a CYP2D6 ultrarapid metabolizer and a CYP2C19 poor metabolizer. This patient may be at an increased risk of a sub-optimal response. Consider selecting an alternative drug.  If Clomipramine is warranted utilize therapeutic drug monitoring to guide dose adjustments. Please consult a clinical pharmacist for more information."
        elif pt[s1] == "normal_metabolizer" and pt[s2] == "ultrarapid_metabolizer":
            recommendation = "This patient is predicted to be a CYP2D6 normal metabolizer and CYP2C19 ultrarapid metabolizer. Based on the CYP2C19 result, the patient may be at an increased risk of a sub-optimal response. Consider selecting an alternative drug not metabolized by CYP2C19. If Clomipramine is warranted utilize therapeutic drug monitoring to guide dose adjustments. Please consult a clinical pharmacist for more information."
        elif pt[s1] == "normal_metabolizer" and pt[s2] == "rapid_metabolizer":
            recommendation = "This patient is predicted to be a CYP2D6 normal metabolizer and CYP2C19 rapid metabolizer. Based on the CYP2C19 result, the patient may be at an increased risk of a sub-optimal response. Consider selecting an alternative drug not metabolized by CYP2C19. If Clomipramine is warranted utilize therapeutic drug monitoring to guide dose adjustments. Please consult a clinical pharmacist for more information."
        elif pt[s1] == "normal_metabolizer" and pt[s2] == "normal_metabolizer":
            recommendation = "No recommendation"
        elif pt[s1] == "normal_metabolizer" and pt[s2] == "intermediate_metabolizer":
            recommendation = "No recommendation"
        elif pt[s1] == "normal_metabolizer" and pt[s2] == "poor_metabolizer":
            recommendation = "This patient is predicted to be a CYP2D6 normal metabolizer  and a CYP2C19 poor metabolizer. Based on the CYP2C19 result, the patient may be at an increased risk of a sub-optimal response. Consider selecting an alternative drug not metabolized by CYP2C19. If Clomipramine is warranted, consider  a 50% reduction of recommended starting dosee and  utilize therapeutic drug monitoring to guide dose adjustments. Please consult a clinical pharmacist for more information."
        elif pt[s1] == "intermediate_metabolizer" and pt[s2] == "ultrarapid_metabolizer":
            recommendation = "This patient is predicated to be a CYP2D6 intermediate metabolizer and a CYP2C19 ultrarapid metabolizer. The patient may be at an increased risk of a sub-optimal response. Consider selecting an alternative drug.  If Clomipramine is warranted utilize therapeutic drug monitoring to guide dose adjustments. Please consult a clinical pharmacist for more information."
        elif pt[s1] == "intermediate_metabolizer" and pt[s2] == "rapid_metabolizer":
            recommendation = "This patient is predicated to be a CYP2D6 intermediate metabolizer and a CYP2C19 rapid metabolizer. The patient may be at an increased risk of a sub-optimal response. Consider selecting an alternative drug.  If Clomipramine is warranted utilize therapeutic drug monitoring to guide dose adjustments. Please consult a clinical pharmacist for more information."
        elif pt[s1] == "intermediate_metabolizer" and pt[s2] == "normal_metabolizer":
            recommendation = "This patient is predicted to be a CYP2D6 intermediate metabolizer and a CYP2C19 normal metabolizer. Based on the CYP2D6 result, the patient may be at an increased risk of a sub-optimal response. Consider  a 25% reduction of recommended starting dosee and  utilize therapeutic drug monitoring to guide dose adjustments."
        elif pt[s1] == "intermediate_metabolizer" and pt[s2] == "intermediate_metabolizer":
            recommendation = "This patient is predicated to be a CYP2D6 intermediate metabolizer and a CYP2C19 intermediate metabolizer. The patient may be at an increased risk of a sub-optimal response. Consider selecting an alternative drug.  If Clomipramine is warranted utilize therapeutic drug monitoring to guide dose adjustments. Please consult a clinical pharmacist for more information."
        elif pt[s1] == "intermediate_metabolizer" and pt[s2] == "poor_metabolizer":
            recommendation = "This patient is predicated to be a CYP2D6 intermediate metabolizer and a CYP2C19 poor metabolizer. The patient may be at an increased risk of a sub-optimal response. Consider selecting an alternative drug.  If Clomipramine is warranted utilize therapeutic drug monitoring to guide dose adjustments. Please consult a clinical pharmacist for more information."
        elif pt[s1] == "poor_metabolizer" and pt[s2] == "ultrarapid_metabolizer":
            recommendation = "This patient is predicated to be a CYP2D6 poor metabolizer and a CYP2C19 ultrarapid metabolizer and may be at an increased risk of a sub-optimal response. Consider selecting an alternative drug.  If Clomipramine is warranted utilize therapeutic drug monitoring to guide dose adjustments. Please consult a clinical pharmacist for more information."
        elif pt[s1] == "poor_metabolizer" and pt[s2] == "rapid_metabolizer":
            recommendation = "This patient is predicated to be a CYP2D6 poor metabolizer and a CYP2C19 rapid metabolizer and may be at an increased risk of a sub-optimal response. Consider selecting an alternative drug.  If Clomipramine is warranted utilize therapeutic drug monitoring to guide dose adjustments. Please consult a clinical pharmacist for more information."
        elif pt[s1] == "poor_metabolizer" and pt[s2] == "normal_metabolizer":
            recommendation = "This patient is predicted to be a CYP2D6  poor metabolizer and a CYP2C19 normal metabolizer. Based on the CYP2D6 result, this patient may be at an increased risk of a sub-optimal response. Consider selecting an alternative drug not metabolized by CYP2D6. If Clomipramine is warranted, consider  a 50% reduction of recommended starting dosee and  utilize therapeutic drug monitoring to guide dose adjustments. Please consult a clinical pharmacist for more information."
        elif pt[s1] == "poor_metabolizer" and pt[s2] == "intermediate_metabolizer":
            recommendation = "This patient is predicated to be a CYP2D6 poor metabolizer and a CYP2C19 intermediate metabolizer . This patient may be at an increased risk of a sub-optimal response. Consider selecting an alternative drug.  If Clomipramine is warranted utilize therapeutic drug monitoring to guide dose adjustments. Please consult a clinical pharmacist for more information."
        elif pt[s1] == "poor_metabolizer" and pt[s2] == "poor_metabolizer":
            recommendation = "This patient is predicated to be a CYP2D6 poor metabolizer and a CYP2C19 poor metabolizer. This patient may be at an increased risk of a sub-optimal response. Consider selecting an alternative drug.  If Clomipramine is warranted utilize therapeutic drug monitoring to guide dose adjustments. Please consult a clinical pharmacist for more information."
        else:
            recommendation = "No recommendation"
    elif s1 in pt and s2 not in pt:
        if pt[s1] == "ultrarapid_metabolizer":
            recommendation = "This patient is predicted to be a CYP2D6 ultrarapid metabolizer and may be at an increased risk of a sub-optimal response. Consider selecting an alternative drug not metabolized by CYP2D6. If Clomipramine is warranted utilize therapeutic drug monitoring to guide dose adjustments. A CYP2C19 genotype does not appear to have been ordered for this patient. CYP2C19 genetic status may be important for alternative drugs. Please consult a clinical pharmacist for more information."
        elif pt[s1] == "normal_metabolizer":
            recommendation = "This patient is predicted to be a CYP2D6 normal metabolizer. There is no reason to selectively adjust the dose of this medication based on the CYP2D6 result. Because a CYP2C19 genotype does not appear to have been ordered for this patient, use of an alternative agent may be recommended. Please consult a clinical pharmacist for more information."
        elif pt[s1] == "intermediate_metabolizer":
            recommendation = "This patient is predicted to be a CYP2D6 intermediate metabolizer and may be at an increased risk of a sub-optimal response. Consider a 25% reduction of recommended starting dosee and  utilize therapeutic drug monitoring to guide dose adjustments. A CYP2C19 genotype does not appear to have been ordered for this patient. Use of an alternative agent may be recommended. Please consult a clinical pharmacistc for more information."
        elif pt[s1] == "poor_metabolizer":
            recommendation = "This patient is predicted to be a CYP2D6 poor metabolizer and may be at an increased risk of a sub-optimal response. Consider selecting an alternative drug not metabolized by CYP2D6. If Clomipramine is warranted, consider  a 50% reduction of recommended starting dosee and  utilize therapeutic drug monitoring to guide dose adjustments. A CYP2C19 genotype does not appear to have been ordered for this patient. CYP2C19 genetic status may be important for alternative drugs. Please consult a clinical pharmacist for more information."
        else:
            recommendation = "No recommendation"
    elif s1 not in pt and s2 in pt:
        if pt[s2] == "ultrarapid_metabolizer":
            recommendation = "This patient is predicted to be a CYP2C19 ultrarapid metabolizer and may be at an increased risk of a sub-optimal response. Consider selecting an alternative drug not metabolized by CYP2C19. If Clomipramine is warranted utilize therapeutic drug monitoring to guide dose adjustments. A CYP2D6 genotype does not appear to have been ordered for this patient. CYP2D6 genetic status may be important for alternative drugs. Please consult a clinical pharmacist for more information."
        if pt[s2] == "rapid_metabolizer":
            recommendation = "This patient is predicted to be a CYP2C19 rapid metabolizer and may be at an increased risk of a sub-optimal response. Consider selecting an alternative drug not metabolized by CYP2C19. If Clomipramine is warranted utilize therapeutic drug monitoring to guide dose adjustments. A CYP2D6 genotype does not appear to have been ordered for this patient. CYP2D6 genetic status may be important for alternative drugs. Please consult a clinical pharmacist for more information."
        if pt[s2] == "normal_metabolizer":
            recommendation = "This patient is predicted to be a CYP2C19 normal metabolizer. There is no reason to selectively adjust the dose of this medication based on the CYP2C19 result. Because a CYP2D6 genotype does not appear to have been ordered for this patient, use of an alternative agent may be recommended. Please consult a clinical pharmacist for more information."
        elif pt[s2] == "intermediate_metabolizer":
            recommendation = "This patient is predicted to be a CYP2C19 intermediate metabolizer. There is no reason to selectively adjust the dose of this medication based on the CYP2C19 result. Because a CYP2D6 genotype does not appear to have been ordered for this patient, use of an alternative agent may be recommended.  Please consult a clinical pharmacist for more information."
        elif pt[s2] == "poor_metabolizer":
            recommendation = "This patient is predicted to be a CYP2C19 poor metabolizer and may be at an increased risk of a sub-optimal response. Consider selecting an alternative drug not metabolized by CYP2C19. If Clomipramine is warranted, consider  a 50% reduction of recommended starting dose and  utilize therapeutic drug monitoring to guide dose adjustments. A CYP2D6 genotype does not appear to have been ordered for this patient. CYP2D6 genetic status may be important for alternative drugs. Please consult a clinical pharmacist for more information."
        else:
            recommendation = "No recommendation"
    elif s1 not in pt and s2 not in pt: 
        recommendation = "No recommendation"
    return recommendation
