# m2_PFA

PFA_fonctionnement.pdf explique comment fonctionne l'algorithme

voici les fichier a ouvrir:

	SDHCAL_PFA.C
	SDHCAL_PFA.h

	Objetcs.C
	Objetcs.h

	Tools.C
	HitGroupingByLayer.C		(ce dernier contient un peu tout c'est la bazar)

(Tools.h et HitGroupingByLayer.h sont vide pour l'instant)

les fichier.root sont dans myroot

Pour lancer un fichier au choix

	$ root *nomDuFichier.root
	root [1] Process("SDHCAL_PFA.C")
	


Le graphe de moi vs Rémi dans PFA_fonctionnement.pdf n'est plus le bon, le bon est results.pdf
sinon pour le produire:

	$ source Macro.sh


