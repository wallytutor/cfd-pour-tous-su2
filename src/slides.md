---
jupyter:
  jupytext:
    cell_metadata_filter: -all
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.18.1
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

# CFD pour tous


## Importation des outils

```python
%load_ext autoreload
%autoreload 2
```

```python
import majordome.latex as ml
```

## Création du *handle*

```python
config = {
    "template": "slides",
    "graphicspath": "../media",
    "title": "CFD pour tous",
    "subtitle": "Introduction à la dynamique des fluides numérique avec SU2",
    "author": "Walter Dal'Maz Silva",
    "date": "\\today"
}

slides = ml.BeamerSlides(config=config, verbose=False)
```

## Introduction

```python
slides.add_section("Introduction")
```

```python
slides.add_slide("intro_prerequis", **{
    "title": "Pré-requis",
    "subtitle": "",
    "contents": ""
})
```

```python
links = {name: ml.url_link(url=url, text=name) for name, url in [
    ("SU2", "https://su2code.github.io/"),
    ("Gmsh", "https://gmsh.info/"),
    ("ParaView", "https://www.paraview.org/"),
    ("Python", "https://www.python.org/"),
    ("PyVista", "https://docs.pyvista.org/"),
    ("MeshLab", "https://www.meshlab.net/"),
    ("Microsoft MPI", "https://learn.microsoft.com/en-us/message-passing-interface/microsoft-mpi"),
]}

with ml.Itemize(itemsep="6pt") as item:
    item.add(f"{links["SU2"]} - le logiciel principal de simulation que nous allons utiliser")
    item.add(f"{links["Gmsh"]} - l'outil capable de générer des mailles pour nos simulations")
    item.add(f"{links["ParaView"]} - le principal outil de visualisation de résults et post-traitement")
    items1 = item.collect()

with ml.Itemize(itemsep="6pt") as item:
    item.add(f"{links["Python"]} - le langage de \\emph{{scripting}} pour l'automatisation de tâches")
    item.add(f"{links["PyVista"]} - la librairie permettant un post-traitement plus avancé des résultats")
    item.add(f"{links["MeshLab"]} - pour une manipulation et correction des maillages génerés")
    item.add(f"{links["Microsoft MPI"]} - pour le calcul parallèle, si votre ordinateur le permet (admin)")
    items2 = item.collect()

with ml.SlideContentWriter() as writer:
    writer.add("Le minimum requis en termes de logiciel :")
    writer.vspace("6mm")
    writer.add(items1)

    writer.vspace("6mm")

    writer.add("Pour les niveaux plus avancés :")
    writer.vspace("6mm")
    writer.add(items2)

    contents = writer.collect()

slides.add_slide("intro_logiciels", **{
    "title": "Logiciels nécessaires",
    "subtitle": "",
    "contents": contents
})
```

## Compilation du document

```python
slides.build("_slides.tex")
```
