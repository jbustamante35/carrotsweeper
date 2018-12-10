study(gaby,roots).
study(nate,roots).
scientist(X):-study(X,roots).
makeSentence-->firstPart,secondPart.
firstPart-->[roots].
firstPart-->[shoots].
secondPart-->[grows].
secondPart-->[dies].
secondPart-->[stops_growing].

flattenOut(_func,_in,_out):-_R=..[_func,_in,_tmp],call(_R),flatten(_tmp,_out).
removeEmpty(_func,_in,_out):-_R=..[_func,_in,_tmp],call(_R),exclude(isEmpty,_tmp,_out).
isEmpty(_in):-_in==[].
geneS(S,GENE):-reverse(S,[d,n,e|TMP]),reverse(TMP,GENE).
searchForGenes([H|T],[_H|_T]):-searchForGenes(T,_T),(peek(H,T,R)->accum(R,[],_H);_H=[]).
searchForGenes([],X):-X=[].
peek(a,[t,g|T],T).
accum([H|T],_inSequence,_gene):-append(_inSequence,[H],_toParse),(geneS(_toParse,_gene)->!;accum(T,_toParse,_gene)).

f(X,Y):-Y is 2*X.
bindX(S,Y):-member(X,S),L=..[f,X,Y],call(L).


genInput(_in,_out,0):-_out=_in.
genInput(_in,_out,_level):-atomic_list_concat([x,_level],_term),(_level > 0->(append([_term],_in,_TMP),_nlevel is _level - 1,genCat(_TMP,_out,_nlevel));[]).

genSingleVar(_level,_out):-atomic_list_concat([x,_level],_out).









