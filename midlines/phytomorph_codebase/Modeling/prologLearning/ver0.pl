%:-use_module(library(file_systems)).
%fail :- 0=1.
hello(world).
animal(dog).
animal(cat).
/*
z(c).
hasa(a,r).
x(A) :- hasa(A,B).
r(B) :- hasa(A,B).


at(a,c).
isR(X,R,Y):-at(X,Y),z(Y),x(X),r(R).

*/
/*
%************************************
% example special case
subClass(X,Y):-directSubClass(X,Y).
subClass(X,Y):-directSubClass(X,Z),subClass(Z,Y).
directSubClass(o1,o2).
directSubClass(o2,o3).

% change notation - think about t(ip,p).
isR(X,ip,Y):-ip(X,Y).
isR(X,p,Y):-isR(X,ip,Y).
isR(X,p,Y):-isR(X,ip,Z),isR(Z,p,Y).
p(X,Y):-isR(X,p,Y).
ip(o1,o2).
ip(o2,o3).


% general attempt with new notation:-fail
isTR(X,W,Y):-isR(X,W,Y),t(W).
isTR(X,W,Y):-isTR(X,IW,Y),t(IW,W).
%isTR(X,w,Y):-isTR(X,iw,Y),t(iw,w).
isTR(X,w,Y):-isTR(X,iw,Z),isTR(Z,w,Y),t(iw,w).
w(X,Y):-isTR(X,w,Y);
t(iw,w).
iw(a1,a2).
iw(a2,a3).
*/
%************************************
%
%************************************
% transitvity:subclassing:domain:range
t(subClass).
type(X,Y):-domain(P,Y),isDR(X,P,_).
type(X,Y):-range(P,Y),isDR(_,P,X).
domain(subClass,class).
range(subClass,class).
%************************************
% make containers:contains:inverse
subClass(logicalContainer,container).
subClass(physicalContainer,container).
t(contains).
t(containedIn).
i(contains,containedIn).
contains(field,box).
contains(box,frogs).
%************************************
% make inheritance of set
isR1(Q,Y):-inher(P),isR(X,P,Y),isR1(Q,X).
isR1(_set,_element):-_R=..[_set,_element],call(_R).
set1(element1).
join(element1,element2).
join(element2,element3).
inher(join).
%************************************
% triple general predicate
% general attempt with general transitive : - true
% implement inverse and transitive inverse
isDR(Y,_iP,X):-i(_P,_iP),isDR(X,_P,Y),!.
isTR(Y,_iP,X):-i(_P,_iP),t(_iP),isTR(X,_P,Y),!.
isDR(X,P,Y):-R =..[P,X,Y],catch(call(R),E,fail).
isTR(X,P,Y):-genTran(X,P,Y),t(P).
genTran(X,P,Y):-isDR(X,P,Y),t(P).
genTran(X,P,Y):-isDR(X,P,Z),genTran(Z,P,Y),t(P).
isR(X,P,Y):-isDR(X,P,Y);isTR(X,P,Y).
%************************************

/*
%************************************
%Test general transitive
t(tw).
t(tq).
w(a1,a2).
tw(o1,o2).
tw(o2,o3).
tq(o1,o2).
tq(o2,o3).
*/
%************************************

%************************************
file('hello').
isDirectory(E):-dirList(S),member(E,S).
%************************************
% example of function apply over list
%**map(func,_input,_output).
map(_func,[_iH|_iT],[_oH|_oT]):-_function=..[_func,_iH,_oH],call(_function),map(_func,_iT,_oT).
map(_,[],[]).
neg(_in,_out):-_out is -_in.
%map(neg,[1,2,3],L).
%************************************
dmap(_func,_constant,[_iH|_iT],[_oH|_oT]):-_function=..[_func,_constant,_iH,_oH],call(_function),dmap(_func,_constant,_iT,_oT).
dmap(_,_,[],[]).
myAppend(_constant,_in,_out):-append([_constant],['/',_in],_tmp),atomic_list_concat(_tmp,_out).
%************************************
abs_Directory_recursive([_iH|_iT],[_oH,_oT]):-directory_directory_abs(_iH,_tmp),abs_Directory_recursive(_tmp,_tmp1),append([_iH],[_tmp1],_oH),abs_Directory_recursive(_iT,_oT).

abs_Directory_recursive([],[]).

directory_directory_abs(_in,_out):-directory_files(_in,_tmp1),dmap(myAppend,_in,_tmp1,_tmp2),include(exists_directory,_tmp2,_tmp3),exclude(endsPeriod,_tmp3,_out).


endsPeriod(String):-removeLastChar(String,_,'.').
dirList(Q):-directory_directory_abs(['/mnt/scratch1/phytomorph_dev'],Q).
%************************************
flattenOut(_func,_in,_out):-_R=..[_func,_in,_tmp],call(_R),flatten(_tmp,_out).
%************************************
assertAsClass([_iH|_iT],Class):-assert(type(_iH,Class)),assertAsClass(_iT,Class).
assertAsClass([],_).
%assertAsClass(hello,type,world).
%************************************
removeLastChar(S1,S2,E):-atom_codes(S1,_tmp),append(S2,[E1],_tmp),char_code(E,E1).
%removeLastChar(S1,Ending,S2):-atomic_list_concat([S2,Ending],S1).
%************************************
% generate and apply function over resulting list : - special of general pattern
genAndmap(_base,_class):-flattenOut(abs_Directory_recursive,_base,_out),write(_out),assertAsClass(_out,_class).
%************************************
name_value(String, Name, Value) :-
        sub_string(String, Before, _, After, '.'), !,
        sub_string(String, 0, Before, _, NameString),
        atom_string(Name, NameString),
        sub_string(String, _, After, 0, ValueString),
	atom_string(Value, ValueString).

%genAndmap(['/mnt/scratch1/phytomorph_dev'],subDir).

/*
%************************************
% example of general reverse notation
%
isR(X,P,Y):-R =..[P,X,Y],call(R).
%************************************
*/
%isR(P,X,Y):-P(X,Y).
%P(X,Y):-isR(X,P,Y).


/*
d(p,S):-p(S,O).
r(p,O):-p(S,O).
map(S,P,O):-d(P,S),r(P,O).
map(S,p,O):-p(S,O);
*/



/*
t(q).
t(w).
w(l1,l2).
w(l2,l3).
pred(w).

*/
/*
p(a,b).
w(a,b).
w(b,c).
t(w).
isR(S,p,O):-p(S,O).

(S,w,O):-w(S,O).
isR(S,w,O):-w(S,O),isR(vO,w,O),t(P).
*/





%itTR(X,p,Y):-p(X,Y),t(p).
%itTR(X,p,Y):-itTR(X,p,Z),itTR(Z,p,Y).




