type t = {
	bait_names: string list;
	bait: GLoc.t;
	other_names: string list;
	other: GLoc.t;
	n_reads: int;
	score: float;
}

let compare c1 c2 =
	match compare c1.bait.chr c2.bait.chr with
		| 0 -> compare c1.bait.lo c2.bait.lo
		| n -> n

let head = "bait_chr\tbait_start\tbait_end\tbait_name\totherEnd_chr\totherEnd_start\totherEnd_end\totherEnd_name\tN_reads\tscore"

let ibed_fmt = ("%s %d %d %s %s %d %d %s %d %f" : _ format)

let ibed_fmt6 = ("%s %d %d %s %s %d %d %s %d %f" : _ format6)

let from_string s =
	Scanf.sscanf s ibed_fmt6 (
		fun b_c b_lo b_hi b_n o_c o_lo o_hi o_n n_read scr -> {
			bait_names = 	String.split_on_char ';' b_n;
			bait = 				GLoc.{chr=b_c; lo=b_lo; hi=b_hi};
			other_names = String.split_on_char ';' o_n;
			other = 			GLoc.{chr=o_c; lo=o_lo; hi=o_hi};
			n_reads = 		n_read;
			score = 			scr;
		}
	)

let to_string c =
	let b_n = String.concat ";" c.bait_names in
	let o_n = String.concat ";" c.other_names in
	Printf.sprintf ibed_fmt
		c.bait.chr	c.bait.lo		c.bait.hi		b_n
		c.other.chr	c.other.lo	c.other.hi	o_n
		c.n_reads	c.score

let sort ?(desc = false) =
	if desc then
		List.sort (fun c1 c2 -> - compare c1 c2)
	else
		List.sort compare

let select_bait (loc : GLoc.t) c_lst =
	let rec aux res = function
		| hd::tl 	when loc.chr <> hd.bait.chr 			-> aux res tl
		| hd::_ 	when loc.hi < hd.bait.lo 					-> res
		| hd::tl 	when GLoc.intersects hd.bait loc 	-> aux (hd::res) tl
		| _::tl 																		-> aux res tl
		| [] 																				-> res
	in
	c_lst
		|> sort
		|> aux []

let select_chr chr c_lst =
	let rec aux res = function
		| hd::tl when res = [] &&
									chr <> hd.bait.chr 	-> aux res tl
		| hd::tl when chr = hd.bait.chr		-> aux (hd::res) tl
		| _ 															-> res
	in
	c_lst
		|> sort
		|> aux []

let group_by_chr c_lst =
	let fold_fun acc c =
		let rest = List.tl acc in
		let prev_grp = List.hd acc in
		let prev_chr = (List.hd prev_grp).bait.chr in
		if c.bait.chr = prev_chr then
			(c::prev_grp)::rest
		else
			[c]::prev_grp::rest
	in
	let c_sorted = sort c_lst in
	List.tl c_sorted
		|> List.fold_left fold_fun [[List.hd c_sorted]]
		|> List.map (fun grp -> ((List.hd grp).bait.chr, grp))

let read_ibed ?(header = true) f_path =
	let idx_offset = if header then 2 else 1 in
	let from_string' idx l =
		try from_string l
		with _ -> failwith (Printf.sprintf "Please check the format of (%s) near line (%d)" f_path (idx + idx_offset))
	in
	let lines = Core.In_channel.read_lines f_path in
	List.mapi from_string' (if header then List.tl lines else lines)

let write_ibed ?(header = true) f_path lst =
	lst
		|> List.map to_string
		|> List.cons (if header then head else "")
		|> Core.Out_channel.write_lines f_path
