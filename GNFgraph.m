%assumption: nodes and edges will be properly initialized
%            nodes are named 1 through n, and inputed in that order

%possible optimizations: preallocate queues and use an index
%                        add ability to handle unequal supply and demand
classdef GNFgraph < handle
    properties
        n
        nodes
        edges
        spanningTree
        flows
    end
    methods
        function obj = GNFgraph()
            obj.n = 0;
            obj.nodes = [];
            obj.edges = [];
            obj.spanningTree = [];
            obj.flows = [];
        end
        function setN(obj,n)
            obj.n = n;
            obj.nodes = NaN*zeros(1,n);
            obj.edges = Inf*ones(obj.n);
            obj.flows = zeros(obj.n);
        end
        function addNode(obj, node, supply)
            %input:
            %   supply - amount of supply (negative if it's a demand node)
            obj.nodes(node) = supply;
        end
        function addEdge(obj, cost, sourceNum, destinationNum)
            %input:
            %   cost - cost to ship 1 unit across the edge
            %   sourceNum - source of edge
            %   destinationNum - destination of edge
            obj.edges(sourceNum,destinationNum) = cost;
        end
        function initSpanningTree(obj)
            %makes initial feasible spanning tree
            %have spanning tree as 2 x (n-1) matrix
            %note: if depleting supply and demand at same time, make sure
            %to add 0 edge
            obj.flows = zeros(obj.n);
            tempNodes = obj.nodes;
            obj.spanningTree = zeros(2,obj.n-1);
            supplys = find(obj.nodes >= 0);
            demands = find(obj.nodes < 0);
            
            %have non-existent edges have very high cost:
            obj.edges(isinf(obj.edges))=1+sum(obj.edges(~isinf(obj.edges)));
            
            %similar to vogel's advanced start method
            %minus the step where if both supply and demand are
            %simultaneously depleted
            for x=1:(obj.n-1)
                %add edge to spanning tree
                
                %find lowest in each column
                if length(supplys) > 1
                    cols = zeros(1,length(demands));
                    for j = 1:length(demands)
                        min2 = Inf;
                        min1 = Inf;
                        for k = supplys
                            if obj.edges(k,demands(j)) < min2
                                min2 = obj.edges(k,demands(j));
                                if min2 < min1
                                    temp = min1;
                                    min1 = min2;
                                    min2 = temp;
                                end
                            end
                        end
                        if isinf(min2)
                            cols(j) = -Inf;
                        else
                            cols(j) = min2 - min1;
                        end
                    end
                else
                    cols = -Inf;
                end
                
                %find lowest in each row
                if length(demands) > 1
                    rows = zeros(1,length(supplys));
                    for j = 1:length(supplys)
                        min2 = Inf;
                        min1 = Inf;
                        for k = demands
                            if obj.edges(supplys(j),k) < min2
                                min2 = obj.edges(supplys(j),k);
                                if min2 < min1
                                    temp = min1;
                                    min1 = min2;
                                    min2 = temp;
                                end
                            end
                        end
                        if isinf(min2)
                            rows(j) = -Inf;
                        else
                            rows(j) = min2 - min1;
                        end
                    end
                else
                    rows = -Inf;
                end
                
                if max(rows) >= max(cols)
                    [temp index] = max(rows);
                    i = supplys(index);
                    min1 = Inf;
                    j = demands(1);
                    for k = demands
                        if obj.edges(i,k) < min1
                            min1 = obj.edges(i,k);
                            j = k;
                        end
                    end
                else
                    [temp index] = max(cols);
                    j = demands(index);
                    min1 = Inf;
                    i = supplys(1);
                    for k = supplys
                        if obj.edges(k,j) < min1
                            min1 = obj.edges(k,j);
                            i = k;
                        end
                    end
                end
                obj.spanningTree(:,x) = [i;j];
                
                obj.flows(i,j) = min(tempNodes(i), -tempNodes(j));
                if tempNodes(i) + tempNodes(j) == 0
                    %remove supply and demand
                    tempNodes(j) = 0;
                    tempNodes(i) = nan;
                    supplys(supplys==i) = [];
                    %demands(demands==j) = [];
                    
                elseif tempNodes(i) + tempNodes(j) > 0
                    %remove demand
                    tempNodes(i) = tempNodes(i) + tempNodes(j);
                    tempNodes(j) = nan;
                    demands(demands==j) = [];
                else
                    %remove supply
                    tempNodes(j) = tempNodes(i) + tempNodes(j);
                    tempNodes(i) = nan;
                    supplys(supplys==i) = [];
                end
            end
        end
        function A = spanningTreeAdjacency(obj)
            A{obj.n,2} = [];
            for i = obj.spanningTree
                A{i(1),1}(end+1) = i(2);
                A{i(2),2}(end+1) = i(1);
            end
        end
        function updateCosts(obj)
            %find y's and have c_i,j = c_i,j - y(i) + y(j)
            y = nan*zeros(1,obj.n);
            adj = spanningTreeAdjacency(obj);
            temp = 1;
            y(temp) = 0;
            queue = [adj{temp,1} adj{temp,2}; temp*ones(1,length(adj{temp,1})...
                +length(adj{temp,2})); ones(1,length(adj{temp,1})) ...
                -ones(1,length(adj{temp,2}))];
            while any(isnan(y))
                temp = queue(:,1);
                queue(:,1) = [];
                if isnan(y(temp(1)))
                    if temp(3) == 1
                        y(temp(1)) = -obj.edges(temp(2),temp(1)) + y(temp(2));
                    else
                        y(temp(1)) = obj.edges(temp(1),temp(2)) + y(temp(2));
                    end
                    temp = temp(1);
                    queue = [queue [adj{temp,1} adj{temp,2}; temp*ones(1,length(adj{temp,1})...
                        +length(adj{temp,2})); ones(1,length(adj{temp,1})) ...
                        -ones(1,length(adj{temp,2}))]];
                end
            end
            obj.edges = obj.edges + ones(obj.n,1)*y - y'*ones(1,obj.n);
        end
        function N = nextEdge(obj)
            %finds next edge to add to the spanning tree
            %(one with negative cost)
            N = [];
            [i j] = find(obj.edges == min(min(obj.edges)),1);
            if obj.edges(i,j) < 0
                N = [i j];
            end
        end
        function shiftFlow(obj, sourceNum, destinationNum)
            %makes specified edge nonbasic
            %changes flow and spanning tree
            %do BFS to get cycle
            BFS{obj.n} = [];
            BFS{destinationNum} = [destinationNum;1];
            adj = spanningTreeAdjacency(obj);
            temp = destinationNum;
            queue = [adj{temp,1} adj{temp,2}; temp*ones(1,length(adj{temp,1})...
                +length(adj{temp,2})); ones(1,length(adj{temp,1})) ...
                -ones(1,length(adj{temp,2}))];
            %make queue: current, previous, direction
            while isempty(BFS{sourceNum})
                temp = queue(:,1);
                queue(:,1) = [];
                if isempty(BFS{temp(1)})
                    curr = temp(1);
                    prev = temp(2);
                    s = temp(3);
                    BFS{curr} = [[curr; s*BFS{prev}(2,1)] BFS{prev}];
                    temp = curr;
                    queue = [queue [adj{temp,1} adj{temp,2}; temp*ones(1,length(adj{temp,1})...
                        +length(adj{temp,2})); s*BFS{prev}(2,1)*[ones(1,length(adj{temp,1})) ...
                        -ones(1,length(adj{temp,2}))]]];
                    %not 100% sure the directions are right
                end
            end
            loop = BFS{sourceNum};
            
            deltaFlow = zeros(4,size(loop,2));
            deltaFlow(:,end) = [sourceNum; destinationNum; 0; 1];
            for i = (size(loop,2)-1):-1:1
                %d = [src; dest; flow; +/-]
                if xor(loop(1,i+1) ~= deltaFlow(1,i+1),...
                        deltaFlow(4,i+1) == loop(2,i))
                    deltaFlow(:,i) = [loop(1,i); loop(1,i+1); ...
                        obj.flows(loop(1,i),loop(1,i+1)); loop(2,i)];
                else
                    deltaFlow(:,i) = [loop(1,i+1); loop(1,i); ...
                        obj.flows(loop(1,i+1),loop(1,i)); loop(2,i)];
                end
            end
            
            cycleFlow = min(deltaFlow(3,deltaFlow(4,:) < 0));
            index = find(deltaFlow(3,:).*deltaFlow(4,:) == -cycleFlow,1);
            
            for i =1:obj.n
                if all(obj.spanningTree(:,i)==deltaFlow(1:2,index))
                    obj.spanningTree(:,i) = [sourceNum; destinationNum];
                    break
                end
            end
            %obj.spanningTree(:,index) = [sourceNum; destinationNum];
            
            for d = deltaFlow
                obj.flows(d(1),d(2)) = obj.flows(d(1),d(2))+cycleFlow*d(4);
            end
        end
        function optimize(obj)
            %optimizes the graph after all input is given
            initSpanningTree(obj);
            while 1
                updateCosts(obj);
                N = nextEdge(obj);
                if isempty(N)
                    break
                end
                shiftFlow(obj, N(1), N(2));
            end
        end
    end
end