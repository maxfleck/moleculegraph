# moleculegraph

simple undirected graph representation for molecules and more...

**Install**

pip install -I .

or

python setup.py install

**moleculegraph**

- package to describe Molecules with focus on force-fields and easy mapping

- main goal: write code about molecules without the need to care abput molecules

- easy semantics, agnositc against syntax

- high flexibillity for naming your beads/atoms/whatever (b and r might have a negative impact. Maybe change)
    
- not limited to small molecules i.e. branches and rings are possible (works for amino acids and therefore proteins too)


**syntax**

<img src="examples/latex/graph_mdma-1.png" width="50%" height="50%">

- [b] encodes the branch operator, [r] encodes the ring operator
- the number after the operator determines its range 
- branches point forward, rings backwards 
- operators refer to the subchain of the last atom before the operator 
- subchains are skipped if completely enclosed by the operator       
- ...so that the number after the operator always corresponds to the number of atoms in the encoded substructure :) 

<img src="examples/latex/graph_mut-1.png" width="50%" height="50%">