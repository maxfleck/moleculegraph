# moleculegraph

simple undirected graph representation for molecules and more...

**Install**

pip install -I .

or

python setup.py install

**moleculegraph**

- class to describe Molecules with focus on force-fields and easy mapping

- main goal: get distance matrices (concept: get them from simple neighbour lists -> works great)

- high flexibillity for Namings (b and r might have a negative impact. Maybe change)
    
- not limited to small molecules i.e. branches and rings are possible


![](latex/mdma.png?raw=true)


**TO DO**

- programming:

    - <span style="color:red"> adjust naming i.e. bond or neighbour ;) </span>
    - <span style="color:red">self.get_bond_matrix etc. is not necessary </span>
    - <span style="color:red">kill self.funs that can crash the class -> ini-stuff is not a function that writes to self!!! </span>
    
- <span style="color:red">merge force field and distance matrix part </span>

- <span style="color:red">make sub-packages i.e. all the defs </span>
