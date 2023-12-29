set Rxn := { read "zimpl_txt/SMatrix_reshaped.txt" as "<2s>" skip 1};
set Met := { read "zimpl_txt/SMatrix_reshaped.txt" as "<3s>" skip 1};
do print card(Met);
set Rev := { read "zimpl_txt/Rev.txt" as "<1s>"};
set Irrev := { read "zimpl_txt/Irrev.txt" as "<1s>"};
set PRxn := { read "zimpl_txt/PRxn.txt" as "<1s>"};
set PMet := { read "zimpl_txt/PMet.txt" as "<1s>"};
set F := { read "zimpl_txt/F.txt" as "<3s>" skip 1};
set Fminus1 := F\{"0"};
set Iterations := { read "zimpl_txt/Iterations.txt" as "<1n>"};
param l[Rxn] :=  read "zimpl_txt/l.txt" as "1n";
param u[Rxn] :=  read "zimpl_txt/u.txt" as "1n";
param delta := read "zimpl_txt/params.txt" as "1n" use 1;
param M := read "zimpl_txt/params.txt" as "1n" skip 1 use 1;
param dof := read "zimpl_txt/params.txt" as "1n" skip 2 use 1;

param S[Rxn*Met] := read "zimpl_txt/SMatrix_reshaped.txt" as "<2s,3s> 4n" skip 1;
param I[Iterations*Rxn] := read "zimpl_txt/prevMinimumSubnetworks.txt" as "<1n,2s> 3n";
param D[Rxn*F] := read "zimpl_txt/F.txt" as "<2s,3s> 4n" skip 1;
param d[F] := read "zimpl_txt/Fd.txt" as "1n";

defset Rxnm(m) := {<i> in Rxn with S[i,m] != 0};
defset Revm(m) := {<i> in Rev with S[i,m] != 0};
defset ActiveRxn(iter) := {<i> in Rxn with I[iter,i] != 0};

var v[Rxn*F] real;
var af[Rxn*F] binary;
var af_bar[Rev*F] binary;
var a[Rxn] binary;
var a_bar[Rev] binary;

# maximize flux: v["Growth", "0"];
minimize network: sum<i> in Rxn: a[i];

subto Functionalities:
      forall <f> in Fminus1:
          sum <i> in Rxn: 
          D[i, f] * v[i, f] <= d[f];

subto SteadyState:
      forall <f> in F:
          forall <i> in Met:
              sum <j> in Rxn:
                  S[j, i] * v[j, f] == 0;

subto Bounds:
      forall <f> in F:
          forall <i> in Rxn:
              l[i] <= v[i, f] and v[i, f] <= u[i];

subto ProtectedRxnIrrev:
      forall <i> in PRxn inter Irrev: 
          af[i, "0"] == 1;

subto ProtectedRxnRev:
      forall <i> in PRxn inter Rev: 
          af[i, "0"] + af_bar[i, "0"] == 1;

subto ProtectedMet:
      forall <m> in PMet do
          if Rxnm(m) inter PRxn == {} then
              if Rev != {} then
                  sum<i> in Rxnm(m): af[i, "0"] + sum<j> in Revm(m): af_bar[j, "0"] >= 1
              else
                  sum<i> in Rxnm(m): af[i, "0"] >= 1
              end
          end;

subto Irrevv:
      forall <f> in F:
          forall <i> in Irrev:
              delta * af[i, f] <= v[i, f] and v[i, f] <= M * af[i, f];

subto Revv:
      forall <f> in F:
          forall <i> in Rev:
              delta * af[i, f]  - M * af_bar[i, f] <= v[i, f] and v[i, f] <= M * af[i, f] - delta * af_bar[i, f];
      
subto Reva:
      forall <f> in F:
          forall <i> in Rev:
              af[i, f] + af_bar[i, f] <= 1;

subto FunIrrev:
      forall <i> in Irrev:
          sum<f> in F: af[i, f] >= a[i] and sum<f> in F: af[i, f] <= a[i] * card(F);

subto FunRev:
      forall <i> in Rev:
          sum<f> in F: (af[i, f] + af_bar[i, f]) >= a[i] and sum<f> in F: (af[i, f] + af_bar[i, f]) <= a[i] * 2 * card(F);

subto prevSubs:
      forall <iter> in Iterations:
          sum <i> in Rxn: 
              (if <i> in ActiveRxn(iter) then (1 - a[i]) else a[i] end) >= 1;
