(**
This program accompanies the paper:
"Cooperation between genetic mutations and phenotypic plasticity can bypass 
the Weismann barrier: The cooperative model of evolution" 
by Ken Nishikawa and Akira R. Kinjo.

This program is written in the OCaml language, version 4.06.1
(http://caml.inria.fr/).

The following third-party softwares are required
in addition to the OCaml compiler.

- Findlib (http://projects.camlcity.org/projects/findlib.html)
- ExtLib (https://code.google.com/p/ocaml-extlib/)
- GSL (https://bitbucket.org/mmottl/gsl-ocaml)

To compile the program, execute the following command:
% ocamlbuild -pkgs bigarray,extlib,gsl adjevol.native

** Running the program **

% ./adjevol -help 
shows available options.

% ./adjevol -seed 313 -theta 5 -qsel 0.15 -ngen 100 -num 100000 -pmut 0.01 -sigma 3 > outfile

reproduces a result of the paper with epigenetic effect (sigma=3).

% ./adjevol -seed 313 -theta 5 -qsel 0.15 -ngen 100 -num 100000 -pmut 0.01 -sigma 0.5 > outfile

reproduces a result of the paper without epigenetic effect. (sigma=0.5).

*)

open ExtLib
open Printf
module Ht = Hashtbl

type individual = {
    genes: int array;
    gef: float;
    epi: float;
    fit: float;
  }

(* epigenetic effect on fitness *)
let epigen setup = 
  let r = setup#epimean +. setup#scale *. setup#sigma *. (Gsl.Randist.ugaussian setup#rng) in
  r 
  (*if r > 0.0 then r else 0.0 *)

(* genetic effect on fitness *)
let fitness_g len w ind =
  let f = ref 0 in
  Array.iter2 (fun w i -> f := !f + w * i) w ind.genes;
  float !f

let make_individual setup w =
  let genes = Array.init setup#len (fun i -> 
    if Gsl.Rng.uniform setup#rng < setup#pmut
    then 1
    else 0) in
  let i = {genes=genes; gef=0.0; epi=0.0; fit=0.0} in
  let g = fitness_g setup#len w i in
  let e = epigen setup in
  {i with gef = g; epi = e; fit = g +. e}

(* set up the weights of genes *)
let make_weights setup = 
  Array.init setup#len (fun i -> 
    if i < setup#len / 2 then -1 else 1)

(* making a population of individuals *)
let make_ensemble setup w =
  Array.init setup#num (fun i -> make_individual setup w)

let setup =
  object
    val rng = Gsl.Rng.make (Gsl.Rng.default ())
    val mutable ngen = 100
    val mutable len = 20
    val mutable theta = 5.0
    val mutable num = 100000
    val mutable max_num = 100000
    val mutable qsel = 0.15
    val mutable pmut = 0.01
    val mutable sigma = 3.0
    val mutable epimean = 0.0
    val mutable scale = 1.0
    method sigma = sigma
    method set_sigma x = sigma <- x
    method scale = scale
    method set_scale x = scale <- x
    method ngen = ngen
    method set_ngen x = ngen <- x
    method rng = rng
    method set_rng i = Gsl.Rng.set rng (Nativeint.of_int i)
    method theta = theta
    method set_theta x = theta <- x
    method len = len
    method set_len x = len <- x
    method num = num
    method set_num x = num <- x
    method max_num = max_num
    method set_max_num x = max_num <- x
    method qsel = qsel
    method pmut = pmut
    method set_qsel x = qsel <- x
    method set_pmut x = pmut <- x
    method epimean = epimean
    method set_epimean x = epimean <- x
  end

(* genetic recombination *)
let mate1 rng len n ens = 
  let i = Gsl.Rng.uniform_int rng n in
  let j = Gsl.Rng.uniform_int rng n in
  let site = Gsl.Rng.uniform_int rng len in
  let ai = {ens.(i) with genes = Array.copy ens.(i).genes} in
  let aj = {ens.(j) with genes = Array.copy ens.(j).genes} in
  for k = 0 to len - 1 do
    if k < site then
      ai.genes.(k) <-ens.(j).genes.(k)
    else
      aj.genes.(k) <- ens.(i).genes.(k)
  done;
  ai,aj

let mate_all rng len ens =
  let n = Array.length ens in
  let n2 = 2*n in
  let nens = Array.create (2*n2) ens.(0) in
  for i = 1 to n2 do
    let a,b = mate1 rng len n ens in
    nens.(2*i - 2) <- a;
    nens.(2*i - 1) <- b
  done;
  nens

let selectp setup ind =
  if ind.fit >= setup#theta then true
  else if ind.fit < 0.0 then false
  else (Gsl.Rng.uniform setup#rng) < setup#qsel

let selectp_prop setup ind =
  if ind.fit >= setup#theta then true
  else if ind.fit < 0.0 then false
  else (Gsl.Rng.uniform setup#rng) < ind.fit /. setup#theta

let selectp_prop2 setup ind =
  let l = 0.5 *. float (setup#len) in
  (Gsl.Rng.uniform setup#rng) < (ind.fit +. l) /. (setup#theta +. l)

let select_all setup w ens =
  let n = Array.length ens - 1 in
  let len = setup#len in
  let lst = ref [] in
  for i = 0 to n do
    let ind = ens.(i) in
    let gcont = fitness_g len w ind in
    let econt = epigen setup in
    let f = gcont +. econt in
    let ind = {ind with gef = gcont; epi = econt; fit=f} in
    if selectp setup ind then lst := ind :: !lst else ()
  done;
  Array.of_list (List.take setup#max_num !lst)

let cycle setup w ens =
  let len = setup#len in
  let rng = setup#rng in
  let ens = mate_all rng len ens in
  select_all setup w ens

let print_individuals setup ens istep = 
  let ht = Ht.create 10 in
  Array.iter (fun i -> Ht.add ht (int_of_float i.gef) ()) ens;
(*
  let nmut = Array.fold_left 
      (fun nmut ind -> 
	if Array.exists ((=) 1) ind.genes then succ nmut else nmut) 
      0 ens in
*)
  let n = Array.length ens in
  fprintf stderr "%d\t%d" istep n;
  (* hoge *)
  let len = setup#len in
  let half_len = len / 2 in
  for i = -half_len to half_len do
    let l = List.length(Ht.find_all ht i) in
    fprintf stderr "\t%d" l
  done;
(*  fprintf stderr "\t%d" nmut;*)
  fprintf stderr "\n";
  flush stderr;
  ()

let print_stats w ens istep = 
  let n = Array.length ens in
  let gtotal = ref 0.0 in
  let etotal = ref 0.0 in
  let ftotal = ref 0.0 in
  Array.iter 
    (fun i -> 
      gtotal := !gtotal +. i.gef;
      etotal := !etotal +. i.epi;
      ftotal := !ftotal +. i.fit)
    ens;
  let fn = float n in
  let gmean = !gtotal /. fn in
  let emean = !etotal /. fn in
  let fmean = !ftotal /. fn in
  let gvar = ref 0.0 in
  let evar = ref 0.0 in
  let fvar = ref 0.0 in
  Array.iter
    (fun i -> 
      let dev = i.gef -. gmean in
      gvar := !gvar +. dev *. dev;
      let dev = i.epi -. emean in
      evar := !evar +. dev *. dev;
      let dev = i.fit -. fmean in
      fvar := !fvar +. dev *. dev)
    ens;
  let gsd = sqrt (!gvar /. fn) in
  let esd = sqrt (!evar /. fn) in
  let fsd = sqrt (!fvar /. fn) in

  printf "step: %5d %5d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n"
    istep n gmean gsd emean esd fmean fsd;
  flush stdout

let print_weights w =
  Array.iteri (fun i w -> printf "# w(%d) = %2d\n" i w) w;
  flush stdout

let print_ens ens =
  let n = Array.length ens in
  let l = Array.length ens.(0).genes in
  for i = 0 to n - 1 do
    print_string "#ens";
    for j = 0 to l - 1 do
      printf "%2d" ens.(i).genes.(j);
    done;
    print_newline()
  done

let _ =
  let specs = [
    ("-seed",Arg.Int (fun i -> setup#set_rng i),"Random seed");
    ("-len",Arg.Int (fun i -> setup#set_len i),"length of chromosome");
    ("-num",Arg.Int (fun i -> setup#set_num i),"population size");
    ("-theta",Arg.Float (fun i -> setup#set_theta i),"threshold value");
    ("-sigma",Arg.Float (fun i -> setup#set_sigma i),"sigma (SD) for epigenetics");
    ("-epimean",Arg.Float (fun i -> setup#set_epimean i),"mean for epigenetics");
    ("-scale",Arg.Float (fun i -> setup#set_scale i),"scaling value for epigenetics");
    ("-ngen",Arg.Int (fun i -> setup#set_ngen i),"Number of generations");
    ("-qsel",Arg.Float (fun i -> setup#set_qsel i),"Probability of random selection");
    ("-pmut",Arg.Float (fun i -> setup#set_pmut i),"Probability of random mutations");
  ] in	       
  Arg.parse specs (fun _ -> ()) "usage:";
  let weights = make_weights setup in
  print_weights weights;

  let ens = make_ensemble setup weights in
(*  
  print_endline "#Initial";
  print_ens ens;
*)

(*  let pw = Array.fold_left (fun p w -> if w > 0 then p +. 1.0 else p) 0.0 weights in*)
(*  setup#set_theta pw;*)
  printf "#threshold = %f\n" setup#theta;
  print_endline "(1) generation, (2) population size, (3) mean genetic effect, (4) S.D. of genetic effect, (5) mean epigenetic effect, (6) S.D. of epigenetic effect, (7) mean fitness, (8) S.D. of fitness";
  print_stats weights ens 0;
  print_individuals setup ens 0;
  let rec loop step ens = 
    if step > setup#ngen then ens
    else
      let ens = cycle setup weights ens in
      if Array.length ens <= 1 then 
	(print_endline "extinct!";
	 exit 0)
      else ();
      print_stats weights ens step;
      print_individuals setup ens step;
      loop (succ step) ens
  in
  let ens = loop 1 ens in
  let _ = ens in
(*
  print_endline "#Final";
  print_ens ens;
*)
  ()
