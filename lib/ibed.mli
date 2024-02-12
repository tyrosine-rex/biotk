type t = {
    bait_names: string list;
    bait: GLoc.t;
    other_names: string list;
    other: GLoc.t;
    n_reads: int;
    score: float;
}

val compare : t -> t -> int

val is_intersect : t -> GLoc.t -> bool

val from_string : string -> t

val to_string : t -> string

val sort : ?desc:bool -> t list -> t list

val select_bait : GLoc.t -> t list -> t list

val select_chr : string -> t list -> t list

val group_by_chr : t list -> (string * t list) list

val read_ibed : ?header:bool -> string -> t list

val write_ibed : ?header:bool -> string -> t list -> unit
