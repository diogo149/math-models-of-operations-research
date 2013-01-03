%GNFsample

GNF = GNFgraph();
GNF.setN(5);
GNF.addNode(1,50);
GNF.addNode(2,-40);
GNF.addNode(3,0);
GNF.addNode(4,10);
GNF.addNode(5,-20);
GNF.addEdge(1,1,3);
GNF.addEdge(5,1,2);
GNF.addEdge(2,2,4);
GNF.addEdge(6,2,5);
GNF.addEdge(3,3,2);
GNF.addEdge(4,3,4);
GNF.addEdge(5,3,5);
GNF.addEdge(2,4,3);
GNF.addEdge(7,4,5);

GNF.optimize