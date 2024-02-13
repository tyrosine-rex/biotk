type item = {
    bait_names: string list;
    bait: GLoc.t;
    other_names: string list;
    other: GLoc.t;
    n_reads: int;
    score: float;
  }

type t = item list

let head = "bait_chr\tbait_start\tbait_end\tbait_name\totherEnd_chr\totherEnd_start\totherEnd_end\totherEnd_name\tN_reads\tscore"

let ibed_fmt = ("%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%d\t%f" : _ format)

let ibed_fmt6 = ("%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%d\t%f" : _ format6)

let orderize it1 it2 =
  GLoc.compare it1.bait it2.bait

let of_line_exn s =
  Scanf.sscanf s ibed_fmt6 (
    fun b_c b_lo b_hi b_n o_c o_lo o_hi o_n n_read scr -> {
      bait_names = String.split_on_char ';' b_n;
      bait = GLoc.{chr=b_c; lo=b_lo; hi=b_hi};
      other_names = String.split_on_char ';' o_n;
      other = GLoc.{chr=o_c; lo=o_lo; hi=o_hi};
      n_reads = n_read;
      score = scr;
    }
  )

let to_line it =
  Printf.sprintf ibed_fmt
    it.bait.chr
    it.bait.lo
    it.bait.hi
    (String.concat ";" it.bait_names)
    it.other.chr
    it.other.lo
    it.other.hi
    (String.concat ";" it.other_names)
    it.n_reads
    it.score

let sort ?(desc = false) =
  List.sort (if not desc then orderize else Fun.flip orderize)

let select_baits_by_loc (loc : GLoc.t) it_lst =
  let rec aux res = function
    | [] -> res
    | hd::tl -> if hd.bait.chr <> loc.chr then (*not same chromosome*)
                  aux res tl
                else if hd.bait.hi < loc.lo then (*same chromosome, but loc is overtake -> break*)
                  res
                else if GLoc.intersects hd.bait loc then (*hd intersect loc*)
                  aux (hd::res) tl
                else (*same chromosome and loc not reached*)
                  aux res tl
  in it_lst
    |> sort ~desc:true
    |> aux []

let select_baits_by_chr chr it_lst =
  let rec aux res = function
    | [] -> res
    | hd::tl -> if res = [] && hd.bait.chr <> chr then
                  aux res tl
                else if hd.bait.chr = chr then
                  aux (hd::res) tl
                else (*not same chromosome but res is filled -> chromosome is overtake -> break*)
                  res
  in it_lst
    |> sort ~desc:true
    |> aux []

let gb_fold_fun comp_func acc (it : item) =
  match acc with
    | [[]] -> [[it]]
    | (hd::tl)::rest -> if comp_func hd it then
                          (it::hd::tl)::rest
                        else
                          [it]::(hd::tl)::rest
    | _ -> failwith "This pattern is not supposed to happen"

let group_by_bait it_lst =
  let fold_fun' = gb_fold_fun (fun c1 c2 -> c1.bait = c2.bait)
  in it_lst
    |> sort
    |> List.fold_left fold_fun' [[]]

let group_by_chr it_lst =
  let fold_fun' = gb_fold_fun (fun c1 c2 -> c1.bait.chr = c2.bait.chr)
  in it_lst
    |> sort
    |> List.fold_left fold_fun' [[]]
    |> List.map (fun grp -> (List.hd grp).bait.chr, grp)

let from_file ?(header = true) f_path =
  let idx_offset = if header then 2 else 1 in
  let from_line' idx l =
    try of_line_exn l
    with _ -> failwith (Printf.sprintf "Please check the format of (%s) near line (%d)" f_path (idx + idx_offset))
  in
  let lines = Core.In_channel.read_lines f_path in
  List.mapi from_line' (if header then List.tl lines else lines)

let to_file ?(header = true) f_path it_lst =
  it_lst
    |> List.map to_line
    |> List.cons (if header then head else "")
    |> Core.Out_channel.write_lines f_path


(* for testing *)


let printf_item it =
  Printf.sprintf "%s --> %s\t (%s) --> (%s)(rd:%d;sc:%f)"
    (GLoc.to_string it.bait)
    (GLoc.to_string it.other)
    (String.concat ";" it.bait_names)
    (String.concat ";" it.other_names)
    it.n_reads it.score

let printf_it_lst it_lst =
  String.concat "\n" @@ List.map printf_item it_lst

let%expect_test "fromfile_sort_orderize_1" =
  Printf.printf "%s" (printf_it_lst @@ sort @@ from_file "../data/test.ibed");
  [%expect {|
    1:100-200 --> 1:800-900	 (A1;A2;A3) --> (B1;B2)(rd:10;sc:1.234500)
    2:2000-3000 --> 2:200-300	 (C1) --> (D1)(rd:100;sc:9.990000)
    2:2300-9300 --> 2:1200-1300	 (C1) --> (D1)(rd:100;sc:9.990000)
    5:1000-1500 --> 5:3000-3500	 (E1) --> (F1)(rd:99;sc:3.141590)
    5:3000-3500 --> 5:60300-60800	 (K1) --> (L1;L2)(rd:729;sc:76.890000)
    5:3200-3890 --> 5:600-800	 (I1) --> (J1;J2)(rd:999;sc:141.590000)
    5:3200-3890 --> 5:70600-70800	 (I1) --> (Q1)(rd:876;sc:19.000000)
    5:3200-4400 --> 5:7300-7800	 (M1;M2) --> (N1)(rd:29;sc:1.590000)
    5:3900-4400 --> 5:11100-11750	 (O1) --> (P1)(rd:29;sc:1.590000)
    X:2000-2500 --> X:12345-71000	 (G1;G2) --> (X1)(rd:3;sc:110.000000)
    X:2000-2500 --> X:13579-81000	 (G1;G2) --> (Y1)(rd:4;sc:120.000000)
    X:2000-2500 --> X:24680-91000	 (G1;G2) --> (Z1)(rd:5;sc:130.000000) |}]

let%expect_test "fromfile_sort_orderize_2" =
  Printf.printf "%s" (printf_it_lst @@ sort ~desc:true @@ from_file "../data/test.ibed");
  [%expect {|
    X:2000-2500 --> X:12345-71000	 (G1;G2) --> (X1)(rd:3;sc:110.000000)
    X:2000-2500 --> X:13579-81000	 (G1;G2) --> (Y1)(rd:4;sc:120.000000)
    X:2000-2500 --> X:24680-91000	 (G1;G2) --> (Z1)(rd:5;sc:130.000000)
    5:3900-4400 --> 5:11100-11750	 (O1) --> (P1)(rd:29;sc:1.590000)
    5:3200-4400 --> 5:7300-7800	 (M1;M2) --> (N1)(rd:29;sc:1.590000)
    5:3200-3890 --> 5:600-800	 (I1) --> (J1;J2)(rd:999;sc:141.590000)
    5:3200-3890 --> 5:70600-70800	 (I1) --> (Q1)(rd:876;sc:19.000000)
    5:3000-3500 --> 5:60300-60800	 (K1) --> (L1;L2)(rd:729;sc:76.890000)
    5:1000-1500 --> 5:3000-3500	 (E1) --> (F1)(rd:99;sc:3.141590)
    2:2300-9300 --> 2:1200-1300	 (C1) --> (D1)(rd:100;sc:9.990000)
    2:2000-3000 --> 2:200-300	 (C1) --> (D1)(rd:100;sc:9.990000)
    1:100-200 --> 1:800-900	 (A1;A2;A3) --> (B1;B2)(rd:10;sc:1.234500) |}]

let%expect_test "fromfile_select_baits_by_loc" =
  let loc = GLoc.({chr="5"; lo=2000; hi=3500}) in
  Printf.printf "%s" (printf_it_lst @@ select_baits_by_loc loc @@ from_file "../data/test.ibed");
  [%expect {|
    5:3000-3500 --> 5:60300-60800	 (K1) --> (L1;L2)(rd:729;sc:76.890000)
    5:3200-3890 --> 5:70600-70800	 (I1) --> (Q1)(rd:876;sc:19.000000)
    5:3200-3890 --> 5:600-800	 (I1) --> (J1;J2)(rd:999;sc:141.590000)
    5:3200-4400 --> 5:7300-7800	 (M1;M2) --> (N1)(rd:29;sc:1.590000) |}]

let%expect_test "fromfile_select_baits_by_chr" =
  Printf.printf "%s" (printf_it_lst @@ select_baits_by_chr "2" @@ from_file "../data/test.ibed");
  [%expect {|
    2:2000-3000 --> 2:200-300	 (C1) --> (D1)(rd:100;sc:9.990000)
    2:2300-9300 --> 2:1200-1300	 (C1) --> (D1)(rd:100;sc:9.990000) |}]

let%expect_test "group_by_bait" =
  let gb = group_by_bait @@ from_file "../data/test.ibed" in
  Printf.printf "%s" (String.concat "\n" @@ List.mapi (fun i c -> Printf.sprintf "\ngrp%d:\n%s" i (printf_it_lst c)) gb);
  [%expect {|
    grp0:
    X:2000-2500 --> X:24680-91000	 (G1;G2) --> (Z1)(rd:5;sc:130.000000)
    X:2000-2500 --> X:13579-81000	 (G1;G2) --> (Y1)(rd:4;sc:120.000000)
    X:2000-2500 --> X:12345-71000	 (G1;G2) --> (X1)(rd:3;sc:110.000000)

    grp1:
    5:3900-4400 --> 5:11100-11750	 (O1) --> (P1)(rd:29;sc:1.590000)

    grp2:
    5:3200-4400 --> 5:7300-7800	 (M1;M2) --> (N1)(rd:29;sc:1.590000)

    grp3:
    5:3200-3890 --> 5:70600-70800	 (I1) --> (Q1)(rd:876;sc:19.000000)
    5:3200-3890 --> 5:600-800	 (I1) --> (J1;J2)(rd:999;sc:141.590000)

    grp4:
    5:3000-3500 --> 5:60300-60800	 (K1) --> (L1;L2)(rd:729;sc:76.890000)

    grp5:
    5:1000-1500 --> 5:3000-3500	 (E1) --> (F1)(rd:99;sc:3.141590)

    grp6:
    2:2300-9300 --> 2:1200-1300	 (C1) --> (D1)(rd:100;sc:9.990000)

    grp7:
    2:2000-3000 --> 2:200-300	 (C1) --> (D1)(rd:100;sc:9.990000)

    grp8:
    1:100-200 --> 1:800-900	 (A1;A2;A3) --> (B1;B2)(rd:10;sc:1.234500) |}]

let%expect_test "group_by_chr" =
  let gb = group_by_chr @@ from_file "../data/test.ibed" in
  Printf.printf "%s" (String.concat "\n" @@ List.map (fun (c, l) -> Printf.sprintf "\n%s:\n%s" c (printf_it_lst l)) gb);
  [%expect {|
    X:
    X:2000-2500 --> X:24680-91000	 (G1;G2) --> (Z1)(rd:5;sc:130.000000)
    X:2000-2500 --> X:13579-81000	 (G1;G2) --> (Y1)(rd:4;sc:120.000000)
    X:2000-2500 --> X:12345-71000	 (G1;G2) --> (X1)(rd:3;sc:110.000000)

    5:
    5:3900-4400 --> 5:11100-11750	 (O1) --> (P1)(rd:29;sc:1.590000)
    5:3200-4400 --> 5:7300-7800	 (M1;M2) --> (N1)(rd:29;sc:1.590000)
    5:3200-3890 --> 5:70600-70800	 (I1) --> (Q1)(rd:876;sc:19.000000)
    5:3200-3890 --> 5:600-800	 (I1) --> (J1;J2)(rd:999;sc:141.590000)
    5:3000-3500 --> 5:60300-60800	 (K1) --> (L1;L2)(rd:729;sc:76.890000)
    5:1000-1500 --> 5:3000-3500	 (E1) --> (F1)(rd:99;sc:3.141590)

    2:
    2:2300-9300 --> 2:1200-1300	 (C1) --> (D1)(rd:100;sc:9.990000)
    2:2000-3000 --> 2:200-300	 (C1) --> (D1)(rd:100;sc:9.990000)

    1:
    1:100-200 --> 1:800-900	 (A1;A2;A3) --> (B1;B2)(rd:10;sc:1.234500) |}]
