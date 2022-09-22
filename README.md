Votre objectif est de créer une librairie de filtres de Bloom capable de charger les k-mers d'un génome puis de le requêter. Les propriétés et requêtes sont expliquées dans la sous partie "Librairie filtre de Bloom".

ATTENTION : Veillez à bien tout lire de ce fichier avant de vous lancer dans la programmation.


# Projet global

Nous vous demandons de produire un logiciel qui prend en ligne de commande:
- 1 fichier au format [FASTA](https://fr.wikipedia.org/wiki/FASTA_(format_de_fichier))
- Un entier k (inférieur ou égal à 31)
- Un entier n pour la taille du filtre de Bloom (max $2^{34}$ (16 Go))
- Un entier nf qui donne le nombre de fonction de hashage (max 64)
- Un entier r pour un nombre de requêtes.

Ce logiciel doit construire en interne un [filtre de Bloom](https://fr.wikipedia.org/wiki/Filtre_de_Bloom) de taille n avec nf fonctions de hashages puis le remplir avec tous les k-mers du fichier passé en entrée.
Il doit ensuite effectuer r requêtes aléatoires is_present sur le filtre de bloom (ie générer r kmers aléatoires et les rechercher)

Exemple de ligne de commande:

```bash
  #         fichier          k    n    nf  r
  ./monprog data/ecoli.fasta 31 456637 3 10000
```

## Attentes

- Un repo Github contenant votre projet dont le lien nous sera envoyé avant Jeudi 6 octobre à midi.
- Un readme dans le Github expliquant comment compiler et exécuter votre code.
- Un code qui compile
- Un code qui s'exécute avec la ligne de commande que vous avez donné.
- Des fichiers commentés
- Un code propre et reprenable par quelqu'un d'autre.
- Langage de programmation: C++ fortement recommandé. Rust ou autre langage bas niveau aussi acceptés.


## Les étapes du projets

Nous vous conseillons de procéder par résolution des étapes suivantes dans l'ordre :  
- Télécharger un génome (voir ci-dessous, commencez par le plus petit)  
- Faites-vous une fonction qui vous permet de lire le fichier fasta de manière à ce que chaque appel à la fonction vous renvoie uniquement le caractère suivant de la séquence ADN  
- Faites une fonction qui encode une chaine de caractères de taille k en une valeur entière (i.e. une fonction de hash).
- Créez un fonction qui à partir d'un kmer précédent et une lecture du caractère suivant dans le fichier vous donne l'entier correspondant au kmer suivant.  
- Créez une classe bloom filter suivant la librairie
- Listez tous les k-mer dans l'ordre du fichier, les hasher, les entrer dans le filtre de bloom
- Codez la fonction de requête du BF et l'utiliser
- Bravo, vous avez fini le code


# Librairie filtre de Bloom

Rappel de ce qu'est un filtre de bloom : https://fr.wikipedia.org/wiki/Filtre_de_Bloom  

Vous devrez coder une librairie (sous la forme de classe) de filtre de Bloom.
Cette classe contient un tableau de Bytes (i.e. 8 bits par Byte).
Si le Bloom filter est de taille n, le tableau fera une taille n/8 + 1.
La classe devra avoir 3 fonctions : un constructeur, `add_value` et `is_present`.

## Constructeur

Le filtre de bloom est construit vide mais avec 2 paramètres:  
- Sa taille (en nombre de bits) n.
- Le nombre de fonctions de hashage utilisées nf.

## add_value

La fonction permet d'ajouter une nouvelle entrée x dans le bloom filter.
Cette entrée doit être hashée nf fois en utilisant la fonction xorshift fournie dans le fichier hash.cpp puis les bits correspondant doivent être mis à 1.
Attention, la valeur 0 ne peut pas être hashée (elle retourne toujours 0).
Trouvez une astuce pour contourner ce problème.

## is_present

La fonction prend en paramètre une valeur x et retourne faux si l'entrée n'est pas dans le BF et vrai si elle l'est.


# Datasets : Séquences ADN

## Téléchargement

Pour télécharger un génome il faut:  
1 - Cliquer sur le lien correspondant au génome (voir ci-après)  
2 - Cliquer en haut à droite sur "send to"  
3 - Sélectionner "File"  
4 - Modifier le format pour mettre "FASTA"  
5 - Cliquer sur "Create File"  

- Virus (Sars-Cov-2): https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2  
- Bactérie (e.coli): https://www.ncbi.nlm.nih.gov/nuccore/NZ_LN874954.1  
- Mamifère (Opossums gris à queue courte): https://www.ncbi.nlm.nih.gov/nuccore/NC_008801.1

## Lecture du fichier .fasta

Les fichiers fasta ci-dessus sont tous formatés de la même manière.
La première ligne est une entête dont le premier caractère est un Chevron fermant.
Les lignes suivantes représentent la séquence ADN.
Les retours à la ligne sont uniquement pour faciliter la lecture d'un Humain et ne représentent rien.
Vous devez donc les ignorer lors de la lecture du fichier.  

Les séquences sont composées principalement des 4 lettres A, C, G, T.
Parfois, les meilleurs machines actuelles n'ont pas réussit à parfaitement lire certains nucléotides.
Ces nucléotides sont remplacés par des N (dans la bactérie et l'humain).
Pour ce projets, ignorez les N.
Ainsi, si la séquence est ACTTNNNNATNGCT considérez à la place que c'est ACTTATGCT.
TTA est donc un 3-mer valide ici.
Votre lecteur de fichier dois passer par dessus ces caractères sans les retourner.
ATTENTION : Lorsque nous testerons votre code, nous utiliserons des séquences avec des N. Votre programme ne dois pas planter.



## k-mers

Pour rappel, un k-mer est un mot de taille k.
Ici, il vous faudra tous les générer à partir de la séquence en entrée.
Ainsi la séquence ACTT avec k=3 doit générer ACT et CTT.
  
Pour pouvoir entrer un k-mer dans un BF, vous devez le transformer en entier.
Pour cela vous pourrez utiliser l'encodage 2 bits de votre choix.
Par exemple avec l'encodage A:0 C:1 T:2 G:3 le kmer CTT vaut $1 * 4^2 + 2 * 4^1 + 2 * 4^0 = 26$.


## Reverse complement.

Comme vous le savez certainement, en génétique, chaque nucléotide a son complément.
Ainsi au sein de la double hélice d'ADN chaque A est pairé avec un T et chaque C avec G.
Ce qui veut alors dire que chaque séquence peut être lue à l'endroit normalement ou à l'envers comme son complément.
ACTT peut donc également être lue comme AAGT.
Chaque k-mer a donc un reverse complement.
Pour choisir entre les deux formes d'un k-mer, nous prendrons toujours celui qui est le plus petit dans l'ordre lexicographique.
Dans l'exemple précédent, la séquence ACTT génère 2 3-mers ACT et CTT.
En appliquant la règle du plus petit lexicographique on génère ACT (plus petit que  son reverse complement AGT) et AAG (plus petit que CTT).
