%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% algorithms: solve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solve( Start, [GH|GP],Solution)  :- reverse(Start,StartTMP),breadthfirst( [StartTMP],GH,GP,[],SolutionTMP),reverse(SolutionTMP,Solution).
%solve([field1],[genoTypeAcc,genoType,v1,b73],S).
breadthfirst( [Path | Paths],GoalStatement,GoalParameters,SolutionsIN,SolutionsFinal):-
	% check if head of potential solutions is a solution
	selectElement(GoalStatement,GoalParameters,[Path],SolutionsIN,SolutionsTMP),
	% if Path is a solution - remove it from itself to be left with empty set
	difference(Path,SolutionsTMP,ToExpand),
	% extent Path (or empty set) to next points that meet traverse goal (hard coded for now)
	extendPath(ToExpand,LongerPaths),
	% filter and accumulate solutions from extended paths
	selectElement(GoalStatement,GoalParameters,LongerPaths,SolutionsTMP,NextLevelSolutions),
	% remove solutions from extended paths to produce residue
	difference(LongerPaths,NextLevelSolutions,LeftOverPotential2),
	% cat residue onto end of Paths
	append(Paths,LeftOverPotential2,NextCall),
	% recurse
	breadthfirst(NextCall, GoalStatement,GoalParameters,NextLevelSolutions,SolutionsFinal).
% base call for breadthfirst. when toSearch is empty, return true.
breadthfirst([],_,_,SolutionsIN,SolutionsIN).
% helper for breadthfirst: extend path to nodes that meet rule contain (hardcoded for now, should change to flex.
extendPath([Node|Path],NewPath):-not(bagof([NewNode,Node|Path],(contain(Node,NewNode),not(member(NewNode,[Node|Path]))),JUNK))->NewPath=[];
	bagof([NewNode,Node|Path],(contain(Node,NewNode),not(member(NewNode,[Node|Path]))),NewPath).
%% apply function over list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n1------n2------n6
% |       |
% |       n9
% |
% n3---------------n5
%      |           |
%      |           |
%      n7          n8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unit test:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% method: insert at position
ins(Val,[H|List],Pos,[H|Res]):- Pos > 1, !, Pos1 is Pos - 1, ins(Val,List,Pos1,Res). 
ins(Val, List, 1, [Val|List]).
% unit test:
%	1:ins(pos2,[pos1,pos3,pos4],2,L).:- [pos1,pos2,pos3,pos4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% method: replace(toReplace,with,inList,Result).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
replace(_, _, [], []).
replace(O, R, [O|T], [R|T2]) :- replace(O, R, T, T2).
replace(O, R, [H|T], [H|T2]) :- H \= O, replace(O, R, T, T2).
% unit test:
%	1: replace(a,b,[a,a],R).:- [b,b].
%	1: replace(*,&,[p,*],R).:- [p,&].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% method: set difference X - Y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
difference(Xs,Ys,D) :- 
  not(findall(X,(member(X,Xs),not(member(X,Ys))),TMP))->D=[];findall(X,(member(X,Xs),not(member(X,Ys))),D).

   



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use: put in set to make transitive :- t(predicate).
% allows transitiveSubPreditcates to assume transitivity.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rule:- isTransitvie if 
% 	1) delcared
% 	2) is transitiveSubPredicateOf
isTransitive(X):-t(X).
isTransitive(X):-transitiveSubPredicateOf(X,Z),isTransitive(Z).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rule:- is direct relationship if
%	1) call is true
%	2) inverse is true via call
isDR(X,P,Y):-_R =..[P,X,Y],catch(call(_R),E,fail).
isDR(X,P,Y):-i(Q,P),isDR(Y,Q,X).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% allow for direct relationship between objects to be
% handled by transitiveSubPredicateOf
% rule:-relationship is true if:
%	1: direct relationship
%	2: inverse relationship
%	3: subPredicate relationship
%	4: transitiveSubPredicate relationship
isSubDR(X,P,Y):-isDR(X,P,Y).
isSubDR(X,P,Y):-i(Q,P),isSubDR(Y,Q,X).
isSubDR(X,P,Y):-subPredicateOf(Q,P),isSubDR(X,Q,Y).
isSubDR(X,P,Y):-transitiveSubPredicateOf(Q,P),isSubDR(X,Q,Y).
% rule:- forward transitive relationship is true
genfTran(X,P,Y):-isSubDR(X,P,Y).
genfTran(X,P,Y):-isSubDR(X,P,Y),isTransitive(P).
genfTran(X,P,Y):-isSubDR(X,P,Z),genfTran(Z,P,Y),isTransitive(P).
% rule:- inverse transitive relationship
geniTran(Y,Q,X):-isSubDR(Y,Q,X).
geniTran(Y,Q,X):-i(P,Q),isSubDR(X,P,Y),isTransitive(P).
geniTran(Y,Q,X):-i(P,Q),isSubDR(X,P,Z),genfTran(Z,P,Y),isTransitive(P).
% rule:- if forward or inverse then transitive is true
genTran(X,P,Y):-genfTran(X,P,Y).
genTran(X,P,Y):-geniTran(X,P,Y).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unit test questions:
%	1: selectElement(longerLAcc,[],[[],[2,2],[2,3,4],[1,2]],[[]],Result).				:-[2,3,4].
%	2: selectElement(longerLAcc,[],[[],[2,2,2],[2,3,4],[1,2]],[[]],Result).				:-[[2, 3, 4], [2, 2, 2]].
% 	3: selectElement(genoTypeAcc,[genoType,v1,b73],[[p1,row1],[p4,what,ever],[p6,no,solution]],[[]],R).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions:
%	1) select object from list: (Goal,Parameters,ListToSearch,currentSelected,ReturnSelectedFromList)
% select element from list
%selectElement(Goal,Parameters,[Head|Tail],Selected) :-
%    selectElement(Goal,Parameters, Tail,[Head],Selected).
% base case of select element from list
selectElement(Goal,Parameters,[],Selected,Selected).
% call statement
selectElement(Goal,Parameters,[Head | Tail], Current, FinalSelected) :-
    call(Goal, Parameters,Head, Current,Selected),
    selectElement(Goal, Parameters,Tail, Selected,FinalSelected).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% accumulators: [PotentialObject,CurrentSolutions,Parameters,ResultList]:given an object,a list of objects, return the one or add to the list and return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	1) get longer list from two lists no extra parameters
longerLAcc(_,L1,[LH|LT],L):-length(L1,LEN1),length(LH,LEN2),(LEN1 > LEN2 -> L = [L1];(LEN1 < LEN2 -> append([LH],LT,L);(append([LH],LT,TMP),append([L1],TMP,L)))).
% 	unit test of longerL:
%	1: longerLAcc([2,2,2],[[1,2,3]],R).			:-[[2, 2, 2], [1, 2, 3]].
% 	2: longerLAcc([2,2,2,2],[[1,2,3]],R).			:-[[2, 2, 2, 2]]
%	3: longerLAcc([2,2],[[1,2,3],[1,2,3]],R).
% 	4: longerLAcc([2,2,2],[[1,2,3],[1,2,3]],R).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	2) cat PotentialObject if meets Function goal with v1 replacement
genoTypeAcc([Function|ParaTail],[PH|PT],Solutions,Result):-replace(v1,PH,ParaTail,ResultParameters),
	apply(Function,ResultParameters) -> 
		append([[PH|PT]],Solutions,Result);
		Result = Solutions.
% cant split [] list therefore define it true and return [].
genoTypeAcc(_,[],_,Result):-Result = [];
% 	unit test of longerL:
%	1: genoTypeAcc([genoType,v1,b73],[p1,row1],[[hello]],R). :-[[p1, row1], [hello]].
%	2: genoTypeAcc([genoType,v1,b73],[[row1],[tail]],[[hello]],R). :-[[p1, row1], [hello]].







predicateChain([],P,F).
predicateChain([H|T],P,F):-predicatePeak(H,T,P,F),predicateChain(T,P,F).
predicatePeak(X,[],_,_).
predicatePeak(X,[H|_],P,F):-_R=..[F,X,P,H],catch(call(_R),_,fail).







t(ancestor).
t(descendant).
i(ancestor,descendant).
i(parent,child).

subPredicateOf(parent,ancestor).
subPredicateOf(maleParent,parent).
subPredicateOf(femaleParent,parent).

maleParent(Father,Child):-cross(_,Father,Child).
femaleParent(Mother,Child):-cross(Mother,_,Child).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make pedigree test
% is1     in1 
% |       | 
% |-fe----ma
% |      |     
% |      f1
% |      |
% |-fe---ma
% |     |
% |     b1
% |     |
% |-fe--ma
%      |
%      b2
% question: given two parents: is1 and in1 --
%		1) have they been crossed
% 		2) how many times been backcrossed - 2
%               3) the list of backcrosses [b1 b2]
%               4) where were b1 and b2 planted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make pedigree test
% p1      p2       p4            	p10
% |       |        |			|
% ---------        |---------------------
%     |            |		   | 	
%     p3           |	           p11		
%     |            |---------------|
%     |-------------               p12
%     |      |     |---------------|	
%     |      p5    |               p13
%     |      |     |---------------|
%     |      |     |               p14
%     |      p6    |
%     |      |-----|
%     --------     |
%            |     |
%            p7    p8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cross(p1,p2,p3).
cross(p3,p4,p5).
cross(p5,p5,p6).
cross(p6,p3,p7).
cross(p6,p4,p8).
cross(p4,p10,p11).
cross(p4,p11,p12).
cross(p4,p12,p13).
cross(p4,p13,p14).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unit test questions:-
%	1: maleParent(p1,p3).		:-true.
%	2: maleParent(p2,p3).		:-false.
%	3: maleParent(p5,p6).		:-true
% 	4: femaleParent(p5,p6). 	:-true
% 	5: genTran(p1,ancestor,p7).	:-true
%	6: genTran(p7,descendant,p1).	:-true
%	7: genTran(p3,descendant,p4). 	:-false
%	8: bagof(X,genTran(p7,ancestor,X),List). :- false
% 	9: setof(X,genTran(X,ancestor,p7),L).	 :- [p1,p2,p3,p4,p5,p6].
% Test parent relationship
% 	10:setof(X,genTran(X,parent,p7),L).	 :- [p3,p6].
% Test child relationship
% 	10:setof(X,genTran(X,child,p4),L).	 :- [p5,p8].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PETIGREE vs FIELDSEASON vs COMBINATION
% where, when,how often reproduced 
% how many times has isolate been backcrossed into inbred


% pattern hunting
backCross(BackCrossParent,[CH|CT]):-isSubDR(BackCrossParent,parent,CH),backCrossAUG(CH,CT),backCross(BackCrossParent,CT).
backCrossAUG(P,[C|_]):-isSubDR(P,parent,C).
backCrossAUG(P,[]).
backCross(BackCrossParent,[]).
backCrossAPI(P1,P2,R):- setof(X,backCross(p4,X),T),selectElement(genoTypeAcc,[isSubDR,P2,parent,v1],T,[],R).
% setof(X,backCross(p4,X),R).
% self pattern

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup for field seasons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% field specific methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
genoTypeInRow(Genotype,Row):-genoType(Genotype),isa(Row,row),contains(Row,Seed),genoType(Seed,Genotype),isa(Seed,seed).




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% declare genotype statements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
genoType(b73).
genoType(mo17).
genoType(cml333).
genoType(ph111).


isa(field,container).
isa(row,container).
isa(plot,container).
isa(station,container).
isa(range,container).

contain(testRow,[s0,s1,s2,s3,s4,s5]).
contain(A,B):-contain(A,L),is_list(L),member(B,L).


harvested(ear,row1).
l(A,B):-[p(A,B),q(A,B)].
l(a,b).
%contain(testRow,s1).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup for field seasons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isa(row1,row).
isa(row2,row).
isa(row3,row).
isa(field1,field).
isa(field2,field).
isa(range1,range).
isa(range2,range).
isa(y2017,year).
isa(y2016,year).


t(contain).
i(contain,containIn).
transitiveSubPredicateOf(genoType,contain).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% declare population genotype(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
genoType(p1,b73).
genoType(p2,mo17).
genoType(p4,b73).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make seeds
isa(p1,seed).
isa(p2,seed).
isa(p4,seed).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make field test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make two fields
contain(y2017,field1).
contain(y2017,field2).
% make two ranges
contain(field1,range1).
contain(field1,range2).
% make three rows
contain(range1,row1).
contain(range1,row2).
contain(range2,row3).
% plant three genotypes
contain(row1,p1).
contain(row2,p2).
contain(row3,p4).


findGenotype(GenoType,InitUnit,FinalUnit):-Goal=[genoType,nV,GenoType],solve(InitUnit,Goal,FinalUnit).


% unit test questions:-
%	1: findGenotype(b73,[field1],Unit).  :-[field1, range1, row1, y2017].





t(hasProperty).


% make containIn the inverse of contain


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make transitiveSubPredicateOf transitive
% make contain transitive
% make genoType a transitiveSubPredicateOfContain
t(transitiveSubPredicateOf).


transitiveSubPredicateOf(containProperty,hasProperty).
transitiveSubPredicateOf(containIn,hasProperty).

contain(s1,e1).
contain(s1,e2).
contain(s2,e3).

hasProperty(s1,pc1).


containProperty(pc1,tall).
containProperty(pc1,purple).
containProperty(pc1,wide).




















