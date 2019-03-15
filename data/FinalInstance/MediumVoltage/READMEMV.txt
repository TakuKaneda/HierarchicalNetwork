This is the medium voltage cluster. Note that at each node where you have a transformer you should connect the low voltage instance I also added.

In this folder I made two files, both are in the csv format:
1) Nodes.txt: In this file each row corresponds to the information of one node of the low voltage network.
      You get respectively (in this exact order) for each node:
          - The integer index of the node
          - It's X position
          - It's Y position
          - 1 if it is a load and 0 otherwise
          - the value in kW of the nominal load (nominal means peak) if the node has an associated load otherwise you get 0
          - 1 if it is the root (of the tree) and 0 otherwise

2) Edges.txt: In this file each row corresponds to one edge of the network:
      You get respectively (in this exact order) for each edge:
          - The integer index of the first node
          - The integer index of the second node

I also added a plot of the sub network, I think it does not need supplementary information. 
