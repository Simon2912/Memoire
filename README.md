# Méthodes numériques - Mémoire

Ce repository contient l’ensemble des codes nécessaires à la simulation de structures en utilisant différentes méthodes numériques.

Le fichier principal à exécuter se trouve dans `src/main.py`. Il permet de lancer la simulation et de choisir la méthode numérique utilisée. Les méthodes disponibles sont : Newton-Raphson par contrôle de la force, la méthode Arc-Length ainsi que la relaxation dynamique implicite. Le choix de la méthode se fait directement dans le fichier `main.py`, tout comme la sélection de la structure à calculer, qui doit être choisie parmi celles présentes dans le dossier `structures_tests`.

Pour effectuer une analyse linéaire, il suffit de définir la variable `linear = True` dans le fichier `main.py`.

Le dossier `DR-hybride` contient l’implémentation de la méthode de relaxation dynamique hybride. Son utilisation consiste à ouvrir le fichier principal du dossier, sélectionner une structure dans `structures_tests`, puis lancer la simulation.

La méthode de relaxation dynamique implémentée dans ce projet est équivalente à celle proposée par le plugin Grasshopper “MUSCLE” disponible ici : https://github.com/JonasFeron/Muscle.git
