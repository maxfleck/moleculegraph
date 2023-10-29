# moleculegraph

[![DOI](https://zenodo.org/badge/673091646.svg)](https://zenodo.org/badge/latestdoi/673091646)

Simple undirected graph representation for molecules and more:

Generates all relevant information about a molecule or the graph of a molecule based on a simple string representation. 
Allows information to be mapped onto the graph to create topologies and other input files, apply force fields to molecular coordinates, determine surrogate model parameters, use group contribution theories and ensures reproducibility.
The aim of this very flexible tool is to be able to write code about molecules without having to worry about the molecules themselves, because moleculegraph takes care of that.


Online documentation is available at [maxfleck.github.io/moleculegraph](https://maxfleck.github.io/moleculegraph/moleculegraph.html)

## installation

pip install -I .

or

python setup.py install

## moleculegraph

- main goal: write code about molecules without the need to care about molecules

- easy syntax, agnositc against semantics

- string representation to describe molecules with focus on force-fields and easy mapping

- high flexibillity for naming your beads/atoms/whatever (b and r might have a negative impact. Maybe change)
    
- not limited to small molecules i.e. branches and rings are possible (works for amino acids and therefore proteins too)


## syntax

<img src="examples/latex/graph_mdma-1.png" width="100%" height="100%">

- [b] encodes the branch operator, [r] encodes the ring operator
- the number after the operator determines its range 
- branches point forward, rings backwards 
- operators refer to the subchain of the last atom before the operator 
- subchains are skipped if completely enclosed by the operator       
- ...so that the number after the operator always corresponds to the number of atoms in the encoded substructure :) 

<img src="examples/latex/graph_mut-1.png" width="100%" height="100%">

## future

- dummy operator to be able to code binding types, for example
- allow different, user definable branch and ring operators
