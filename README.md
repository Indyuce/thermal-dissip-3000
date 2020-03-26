# thermal-dissip-3000
Modelisation de transfert thermique par conduction dans un milieu isotrope et uniforme

Lien utile: https://fr.wikipedia.org/wiki/Conduction_thermique

Conductivité thermique massique de l'aluminium à 20°C: 185 W/mK\
Notation: lambda

On modélise le solide métallique par un quadrillage de points entre lesquels il y a une résistance thermique.\
https://fr.wikipedia.org/wiki/R%C3%A9sistance_thermique_de_conduction \
Flux de chaleur: en W i.e J/s, notation: PHI majuscule

PHI = dT / R

La résistance thermique s'exprime en fonction de l'épaisseur de la résistance, de la surface de la résistance, et de la conductivité thermique. La résistance est inversement proportionelle à la conductivité et à la surface, proportionelle à l'epaisseur.

R = e / (lambda * S)

Nécessité de passer par l'énergie thermique pour calculer un transfert de température de proche en proche. Lien température-énergie thermique: https://fr.wikipedia.org/wiki/%C3%89nergie_thermique

Q = m * c * dT

Q = transfert d'énergie thermique (différence d'énergie, en J)\
m = masse du solide considéré\
c = capacité thermique massique du solide

On a aussi l'expression du flux de chaleur:

phi = Q / dt\
Donc Q = phi * dt\
Q = (dT / R) * dt\
Q = dT / (e / lambda * S) * dt\
Q = dT * (dt * lambda * S / e)

Maintenant on cherche la différence de température transmise.

Q = mc * dTt
Donc dTt = dT * (dt * lambda * S / e / mc)

(dt * lambda * S / e / mc) est une constante donc on la calcule avant pour limiter les calculs à chaque itération. Plus dt est petit plus cette approximation est bonne.
Sorte de méthode d'Euler pour avoir une solution approchée.

La température transmise de proche en proche dTt pendant dt est proportionelle à la différence de température actuelle dT, avec un coefficient de proportionnalité constant qui dépend des
caractéristiques géométriques du dissipateur, du matériau utilisé...
