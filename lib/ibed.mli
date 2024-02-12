type t = {
    bait_names: string list;
    bait: GLoc.t;
    other_names: string list;
    other: GLoc.t;
    n_reads: int;
    score: float;
}

val compare : t -> t -> int

val from_line : string -> t

val to_line : t -> string

val sort : ?desc:bool -> t list -> t list

val select_baits_by_loc : GLoc.t -> t list -> t list

val select_baits_by_chr : string -> t list -> t list

val group_by_bait : t list -> t list list

val group_by_chr : t list -> (string * t list) list

val read_ibed : ?header:bool -> string -> t list

val write_ibed : ?header:bool -> string -> t list -> unit
