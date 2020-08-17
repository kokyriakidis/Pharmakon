def Sinvastatin_recommendation(pt: str) -> str:
    recommendation = "slco1b1"
    s1 = ""
    if s1 in pt:
        if pt[s1] == "normal_function":
            recommendation = "Prescribe desired starting dose and adjust doses of Simvastatin based on disease-specific guidelines.	"
        elif pt[s1] == "decreased_function":
            recommendation = "Prescribe a lower dose or consider an alternative statin (e.g. Pravastatin or Rosuvastatin); consider routine CK surveillance."
        elif pt[s1] == "poor_function":
            recommendation = "Prescribe a lower dose or consider an alternative statin (e.g. Pravastatin or Rosuvastatin); consider routine CK surveillance."
        else:
            recommendation = "No recommendation"
    elif s1 not in pt: 
        recommendation = "No recommendation"
    return recommendation