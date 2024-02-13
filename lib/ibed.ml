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
