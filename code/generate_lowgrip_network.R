# ═══════════════════════════════════════════════════════════════════
# Low Grip Strength integrated 3D multimorbidity network
#
# 保留原四层设计：
# Center + Upstream positive + Downstream positive
#        + Up-associated + Down-associated
#
# 新增：
# Cross-associated（中心前后较近）
# Background Community network（中心前后较远）
#
# Disease universe:
# 283 Community diseases
# + M45/L60/K74等network外单向trajectory-positive diseases
# = 286 diseases
# + Low Grip Strength center
# = 287 nodes
#
# Edge hierarchy:
# trajectory_up/down       = temporal LGS trajectory
# association_up/down      = differential edge around trajectory disease
# association_cross        = differential edge to cross-associated disease
# network_background       = underlying LowGrip multimorbidity network
# ═══════════════════════════════════════════════════════════════════

rm(list=ls());gc()
library(tidyverse);library(data.table);library(jsonlite);library(igraph)

# ── 0. 路径 ───────────────────────────────────────────────
direction_path <- "F:/0-共病分析-UKBIOBANK/UKBIOBANK变量数据/Raw.data/Final_Grip_strength_cohort/Results结果/Part3-轨迹分析/Direction_Results/"

# ⭐ 修改成你实际Network_Comparison目录
network_compare_path <- "F:/0-共病分析-UKBIOBANK/UKBIOBANK变量数据/Raw.data/Final_Grip_strength_cohort/Results结果/Part2-网络构建-网络对比/Part1-网络对比-Network_Results_lowgrip_20260706/Network_Comparison/"

output_path <- paste0(direction_path,"JSON_Output_LowGrip_4layer_FullNetwork/")
dir.create(output_path,showWarnings=FALSE,recursive=TRUE)

pick_file <- function(paths,label){
  hit <- paths[file.exists(paths)]
  if(length(hit)==0) stop("❌ 未找到",label,"\n候选路径:\n",paste(paths,collapse="\n"))
  hit[1]
}

community_file <- pick_file(c(
  paste0(network_compare_path,"Community_Membership_res1.0.csv"),
  paste0(network_compare_path,"Community_Membership.csv"),
  paste0(direction_path,"Community_Membership_res1.0.csv"),
  paste0(direction_path,"Community_Membership.csv")
),"Community Membership")

edge_file <- pick_file(c(
  paste0(network_compare_path,"Edge_Comparison_Summary_FDR.csv"),
  paste0(direction_path,"Edge_Comparison_Summary_FDR.csv")
),"Edge Comparison")

lowgrip_network_file <- pick_file(c(
  paste0(network_compare_path,"LowGrip_Common_Network.rds"),
  paste0(direction_path,"LowGrip_Common_Network.rds")
),"LowGrip Common Network RDS")

mapping_file <- pick_file(c(
  paste0(direction_path,"combinedICD.csv"),
  "F:/0-共病分析-UKBIOBANK/UKBIOBANK变量数据/Raw.data/X_Zdata/combinedICD.csv"
),"combinedICD")

engine_A_file <- paste0(direction_path,"EngineA_Upstream_Clogit_LowGrip.csv")
engine_B_file <- paste0(direction_path,"EngineB_Downstream_Cox_LowGrip.csv")

if(!file.exists(engine_A_file)) stop("❌ Engine A文件不存在")
if(!file.exists(engine_B_file)) stop("❌ Engine B文件不存在")

# ── 1. 辅助函数与配色 ─────────────────────────────────────
clean_code <- function(x){
  x <- toupper(trimws(as.character(x)))
  x[x==""] <- NA_character_
  x
}

first_nonempty <- function(x){
  x <- as.character(x)
  x <- x[!is.na(x)&x!=""]
  if(length(x)==0) NA_character_ else x[1]
}

category_color_map <- c(
  "Central"="#1F618D",
  "Infectious disease"="#FFA01A",
  "Circulatory system disease"="#FF8C42",
  "Dermatologic disease"="#FFE66D",
  "Digestive system disease"="#A8E6CF",
  "Diseases of sense organs"="#87CEEB",
  "Diseases of the nervous system"="#B0C4DE",
  "Endocrine/metabolic disease"="#FFB6C1",
  "Genitourinary system disease"="#DDA0DD",
  "Hematopoietic disease"="#98D8C8",
  "Mental disorders"="#C8A2C8",
  "Musculoskeletal system disease"="#F7DC6F",
  "Respiratory system disease"="#90EE90",
  "Neoplasms"="#E74C3C",
  "Sarcopenia"="#2E86C1",
  "Other/Unclassified"="#CCCCCC"
)

get_color <- function(x){
  out <- unname(category_color_map[match(x,names(category_color_map))])
  out[is.na(out)] <- "#CCCCCC"
  out
}

pair_key <- function(a,b){
  a <- as.character(a);b <- as.character(b)
  paste(pmin(a,b),pmax(a,b),sep="__")
}

# ── 2. ICD分类映射 ────────────────────────────────────────
combined_mapping <- fread(mapping_file) %>%
  select(
    Disease_Code=`Combined ICD-10 Codes`,
    `Disease category`
  ) %>%
  mutate(Disease_Code=clean_code(Disease_Code)) %>%
  filter(!is.na(Disease_Code)) %>%
  distinct(Disease_Code,.keep_all=TRUE)

# ── 3. Community Membership：正式283-node network universe ─
community <- read_csv(community_file,show_col_types=FALSE)

if(!"Disease_Name" %in% names(community))
  community$Disease_Name <- NA_character_
if(!"Category" %in% names(community))
  community$Category <- NA_character_
if(!"LowGrip_Module" %in% names(community))
  community$LowGrip_Module <- NA_character_

community <- community %>%
  mutate(
    Disease_Code=clean_code(Disease_Code),
    Disease_Name=as.character(Disease_Name),
    Category=as.character(Category),
    LowGrip_Module=as.character(LowGrip_Module)
  ) %>%
  filter(!is.na(Disease_Code)) %>%
  distinct(Disease_Code,.keep_all=TRUE) %>%
  left_join(combined_mapping,by="Disease_Code") %>%
  mutate(
    Disease_Name=coalesce(Disease_Name,Disease_Code),
    Category=coalesce(Category,`Disease category`,"Other/Unclassified")
  ) %>%
  select(Disease_Code,Disease_Name,Category,LowGrip_Module,everything())

community_codes <- unique(community$Disease_Code)

cat("\n====== COMMUNITY MEMBERSHIP ======\n")
cat("Community diseases:",length(community_codes),"\n")
cat("Duplicate codes:",sum(duplicated(community_codes)),"\n")

if(length(community_codes)!=283)
  stop("❌ Community疾病数不是283")

# ── 4. 完整LowGrip common network ─────────────────────────
lowgrip_net <- readRDS(lowgrip_network_file)
V(lowgrip_net)$name <- clean_code(V(lowgrip_net)$name)
network_rds_codes <- unique(V(lowgrip_net)$name)

cat("\n====== LOWGRIP COMMON NETWORK ======\n")
cat("RDS nodes:",length(network_rds_codes),"\n")
cat("RDS edges:",ecount(lowgrip_net),"\n")
cat("Community但RDS没有:",length(setdiff(community_codes,network_rds_codes)),"\n")
cat("RDS但Community没有:",length(setdiff(network_rds_codes,community_codes)),"\n")

if(!setequal(community_codes,network_rds_codes))
  stop("❌ Community Membership与LowGrip_Common_Network节点不一致")

# ── 5. Engine A/B：Bonferroni真阳性 ───────────────────────
engine_A_all <- fread(engine_A_file) %>%
  rename(Disease_Code=disease) %>%
  mutate(Disease_Code=clean_code(Disease_Code),engine="A_Upstream") %>%
  distinct(Disease_Code,.keep_all=TRUE)

engine_B_all <- fread(engine_B_file) %>%
  rename(Disease_Code=disease) %>%
  mutate(Disease_Code=clean_code(Disease_Code),engine="B_Downstream") %>%
  distinct(Disease_Code,.keep_all=TRUE)

req <- c("Disease_Code","estimate","p_bonferroni")
if(length(setdiff(req,names(engine_A_all)))>0) stop("❌ Engine A缺必要字段")
if(length(setdiff(req,names(engine_B_all)))>0) stop("❌ Engine B缺必要字段")

engine_A <- engine_A_all %>%
  filter(
    !is.na(p_bonferroni),
    p_bonferroni<0.05,
    !is.na(estimate),
    is.finite(estimate),
    estimate>1
  )

engine_B <- engine_B_all %>%
  filter(
    !is.na(p_bonferroni),
    p_bonferroni<0.05,
    !is.na(estimate),
    is.finite(estimate),
    estimate>1
  )

A_positive_codes <- unique(engine_A$Disease_Code)
B_positive_codes <- unique(engine_B$Disease_Code)

bidirectional <- intersect(A_positive_codes,B_positive_codes)
up_codes <- setdiff(A_positive_codes,B_positive_codes)
down_codes <- setdiff(B_positive_codes,A_positive_codes)
all_positive_codes <- unique(c(up_codes,down_codes))

cat("\n====== TRAJECTORY CLASSIFICATION ======\n")
cat("Engine A true positive:",length(A_positive_codes),"\n")
cat("Engine B true positive:",length(B_positive_codes),"\n")
cat("Upstream only:",length(up_codes),"\n")
cat("Downstream only:",length(down_codes),"\n")
cat("Bidirectional:",length(bidirectional),"\n")

if(length(up_codes)!=25) warning("⚠️ Upstream不是预期25")
if(length(down_codes)!=41) warning("⚠️ Downstream不是预期41")
if(length(bidirectional)!=10) warning("⚠️ Bidirectional不是预期10")

# ── 6. Community vs trajectory ────────────────────────────
up_in_network <- intersect(up_codes,community_codes)
down_in_network <- intersect(down_codes,community_codes)
bi_in_network <- intersect(bidirectional,community_codes)

up_outside <- setdiff(up_codes,community_codes)
down_outside <- setdiff(down_codes,community_codes)
bi_outside <- setdiff(bidirectional,community_codes)

trajectory_only <- unique(c(up_outside,down_outside))

visual_disease_codes <- unique(c(
  community_codes,
  up_codes,
  down_codes
))

cat("\n====== NETWORK vs TRAJECTORY ======\n")
cat("Community:",length(community_codes),"\n")
cat("Upstream in Community:",length(up_in_network),"\n")
cat("Upstream outside:",length(up_outside)," -> ",paste(up_outside,collapse=", "),"\n")
cat("Downstream in Community:",length(down_in_network),"\n")
cat("Downstream outside:",length(down_outside)," -> ",paste(down_outside,collapse=", "),"\n")
cat("Bidirectional in Community:",length(bi_in_network),"\n")
cat("Bidirectional outside:",length(bi_outside)," -> ",paste(bi_outside,collapse=", "),"\n")
cat("Network外单向trajectory:",length(trajectory_only)," -> ",paste(trajectory_only,collapse=", "),"\n")
cat("最终疾病universe:",length(visual_disease_codes),"\n")

if(length(visual_disease_codes)!=286)
  warning("⚠️ 最终疾病universe不是预期286")

# ── 7. 名称映射 ───────────────────────────────────────────
get_engine_names <- function(x){
  if("Disease_Name" %in% names(x)){
    x %>% transmute(Disease_Code,Disease_Name=as.character(Disease_Name))
  }else{
    tibble(Disease_Code=x$Disease_Code,Disease_Name=NA_character_)
  }
}

engine_name_map <- bind_rows(
  get_engine_names(engine_A_all),
  get_engine_names(engine_B_all)
) %>%
  group_by(Disease_Code) %>%
  summarise(Disease_Name_engine=first_nonempty(Disease_Name),.groups="drop")

centrality_candidates <- c(
  paste0(network_compare_path,"Centrality_Comparison_Hub.csv"),
  paste0(direction_path,"Centrality_Comparison_Hub.csv")
)

centrality_file <- centrality_candidates[file.exists(centrality_candidates)][1]

if(length(centrality_file)>0&&!is.na(centrality_file)){
  node_name_map <- read_csv(centrality_file,show_col_types=FALSE) %>%
    transmute(
      Disease_Code=clean_code(Disease_Code),
      Disease_Name_centrality=as.character(Disease_Name)
    ) %>%
    distinct(Disease_Code,.keep_all=TRUE)
}else{
  node_name_map <- tibble(
    Disease_Code=character(),
    Disease_Name_centrality=character()
  )
}

make_meta <- function(codes){
  tibble(Disease_Code=unique(codes)) %>%
    left_join(
      community %>%
        select(
          Disease_Code,
          Disease_Name_community=Disease_Name,
          Category_community=Category,
          LowGrip_Module
        ),
      by="Disease_Code"
    ) %>%
    left_join(node_name_map,by="Disease_Code") %>%
    left_join(engine_name_map,by="Disease_Code") %>%
    left_join(combined_mapping,by="Disease_Code") %>%
    mutate(
      Disease_Name=coalesce(
        Disease_Name_community,
        Disease_Name_centrality,
        Disease_Name_engine,
        Disease_Code
      ),
      Category=coalesce(
        Category_community,
        `Disease category`,
        "Other/Unclassified"
      ),
      LowGrip_Module=as.character(LowGrip_Module)
    ) %>%
    select(Disease_Code,Disease_Name,Category,LowGrip_Module)
}

up_traj <- engine_A %>%
  filter(Disease_Code %in% up_codes) %>%
  select(Disease_Code,estimate) %>%
  distinct(Disease_Code,.keep_all=TRUE) %>%
  left_join(make_meta(up_codes),by="Disease_Code")

down_traj <- engine_B %>%
  filter(Disease_Code %in% down_codes) %>%
  select(Disease_Code,estimate) %>%
  distinct(Disease_Code,.keep_all=TRUE) %>%
  left_join(make_meta(down_codes),by="Disease_Code")

# ── 8. Differential edge数据 ──────────────────────────────
edge_data <- read_csv(edge_file,show_col_types=FALSE) %>%
  mutate(
    Disease1_Code=clean_code(Disease1_Code),
    Disease2_Code=clean_code(Disease2_Code)
  )

edge_codes <- unique(c(
  edge_data$Disease1_Code,
  edge_data$Disease2_Code
))
edge_codes <- edge_codes[!is.na(edge_codes)]

cat("\n====== EDGE COMPARISON COVERAGE ======\n")
cat("Edge file diseases:",length(edge_codes),"\n")
cat("Community但Edge没有:",length(setdiff(community_codes,edge_codes)),"\n")

fdr_col <- grep(
  "fdr|q_value",
  names(edge_data),
  value=TRUE,
  ignore.case=TRUE
)[1]

if(is.na(fdr_col))
  stop("❌ Edge Comparison文件没有FDR/q-value列")

# ── 9. 原来的第四层关联疾病 ───────────────────────────────
# Bidirectional不允许进入左右associated层，
# 它们稍后作为Community background展示
excluded_assoc_targets <- unique(c(
  all_positive_codes,
  bidirectional,
  "SARCOPENIA",
  "LOW_GRIP_STRENGTH"
))

build_fourth <- function(pos_codes){
  edges <- edge_data %>%
    filter(
      !is.na(.data[[fdr_col]]),
      .data[[fdr_col]]<0.05
    ) %>%
    filter(
      (Disease1_Code %in% pos_codes &
         !Disease2_Code %in% excluded_assoc_targets) |
        (Disease2_Code %in% pos_codes &
           !Disease1_Code %in% excluded_assoc_targets)
    ) %>%
    mutate(
      from=ifelse(
        Disease1_Code %in% pos_codes,
        Disease1_Code,
        Disease2_Code
      ),
      to=ifelse(
        Disease1_Code %in% pos_codes,
        Disease2_Code,
        Disease1_Code
      )
    ) %>%
    filter(!to %in% excluded_assoc_targets) %>%
    distinct(from,to)
  
  others <- unique(edges$to)
  
  nodes <- make_meta(others) %>%
    transmute(
      name=Disease_Code,
      display_name=Disease_Name,
      category=Category,
      LowGrip_Module
    )
  
  list(edges=edges,nodes=nodes)
}

fourth_up_all <- build_fourth(up_codes)
fourth_down_all <- build_fourth(down_codes)

# 同时与上下游阳性疾病存在differential edge的节点
cross_nodes <- intersect(
  fourth_up_all$edges$to,
  fourth_down_all$edges$to
)

# ⭐ 保持原四层设计：
# cross从左右associated层移走，但不删除其真实边
fourth_up <- list(
  edges=fourth_up_all$edges %>% filter(!to %in% cross_nodes),
  nodes=fourth_up_all$nodes %>% filter(!name %in% cross_nodes)
)

fourth_down <- list(
  edges=fourth_down_all$edges %>% filter(!to %in% cross_nodes),
  nodes=fourth_down_all$nodes %>% filter(!name %in% cross_nodes)
)

cat("\n====== ORIGINAL FOUR-LAYER STRUCTURE ======\n")
cat("Up-associated:",nrow(fourth_up$nodes),"\n")
cat("Down-associated:",nrow(fourth_down$nodes),"\n")
cat("Cross-associated:",length(cross_nodes),"\n")

# ── 10. 原四层坐标：保持原设计与随机调用顺序 ──────────────
set.seed(123)

x_pos <- 17.5
y_range_pos <- 50
x_assoc_min <- 30
x_assoc_max <- 40
y_range_assoc <- 55
z_range <- 15

make_z_seq <- function(n,r=z_range){
  if(n<=0) return(numeric(0))
  if(n==1) return(sample(c(-1,1),1)*runif(1,1,5))
  z <- seq(-r,r,length.out=n)
  if(0 %in% round(z,6))
    z <- z+diff(z[1:2])/2
  z
}

center_id <- "LOW_GRIP_STRENGTH"

center_node <- tibble(
  id=center_id,
  label="Low Grip Strength",
  category="Central",
  color="#1F618D",
  module=NA_character_,
  network_member=FALSE,
  x=0,y=0,z=0,
  size=1.0,
  group="center"
)

# 原上游阳性：完全保持左侧
up_pos <- up_traj %>%
  mutate(
    id=Disease_Code,
    label=Disease_Name,
    category=Category,
    color=get_color(category),
    module=LowGrip_Module,
    network_member=id %in% community_codes,
    x=-x_pos,
    y=seq(y_range_pos,-y_range_pos,length.out=n()),
    z=0,
    size=1.0,
    group="up_positive"
  ) %>%
  select(
    id,label,category,color,module,network_member,
    x,y,z,size,group
  )

# 原下游阳性：完全保持右侧
down_pos <- down_traj %>%
  mutate(
    id=Disease_Code,
    label=Disease_Name,
    category=Category,
    color=get_color(category),
    module=LowGrip_Module,
    network_member=id %in% community_codes,
    x=x_pos,
    y=seq(y_range_pos,-y_range_pos,length.out=n()),
    z=0,
    size=1.0,
    group="down_positive"
  ) %>%
  select(
    id,label,category,color,module,network_member,
    x,y,z,size,group
  )

# 原up-associated：保持原随机坐标逻辑
n_up <- nrow(fourth_up$nodes)

up_assoc <- if(n_up>0){
  fourth_up$nodes %>%
    mutate(
      id=name,
      label=display_name,
      category=coalesce(category,"Other/Unclassified"),
      color=get_color(category),
      module=LowGrip_Module,
      network_member=TRUE,
      x=runif(n(),-x_assoc_max,-x_assoc_min),
      y=runif(n(),-y_range_assoc,y_range_assoc),
      z=make_z_seq(n_up),
      size=1.0,
      group="up_associated"
    ) %>%
    select(
      id,label,category,color,module,network_member,
      x,y,z,size,group
    )
}else tibble()

# 原down-associated：保持原随机坐标逻辑
n_down <- nrow(fourth_down$nodes)

down_assoc <- if(n_down>0){
  fourth_down$nodes %>%
    mutate(
      id=name,
      label=display_name,
      category=coalesce(category,"Other/Unclassified"),
      color=get_color(category),
      module=LowGrip_Module,
      network_member=TRUE,
      x=runif(n(),x_assoc_min,x_assoc_max),
      y=runif(n(),-y_range_assoc,y_range_assoc),
      z=make_z_seq(n_down),
      size=1.0,
      group="down_associated"
    ) %>%
    select(
      id,label,category,color,module,network_member,
      x,y,z,size,group
    )
}else tibble()

# ── 11. 新增：LowGrip中心前后网络层 ────────────────────────
original_front_codes <- unique(c(
  up_pos$id,
  down_pos$id,
  up_assoc$id,
  down_assoc$id
))

# cross来自common network，应属于Community
cross_codes <- intersect(
  cross_nodes,
  community_codes
)

cross_outside <- setdiff(
  cross_nodes,
  community_codes
)

if(length(cross_outside)>0){
  cat("⚠️ Cross但不属于Community:\n")
  print(cross_outside)
}

# Community中尚未被原四层使用的全部疾病
# 包括9个Community内bidirectional
background_codes <- setdiff(
  community_codes,
  unique(c(original_front_codes,cross_codes))
)

cat("\n====== NEW Z-AXIS NETWORK CONTEXT ======\n")
cat("原四层疾病:",length(original_front_codes),"\n")
cat("Cross-associated:",length(cross_codes),"\n")
cat("Background Community:",length(background_codes),"\n")

# cross：靠近中心前后
cross_meta <- make_meta(cross_codes)
n_cross <- nrow(cross_meta)

cross_nodes_3d <- if(n_cross>0){
  cross_meta %>%
    mutate(
      id=Disease_Code,
      label=Disease_Name,
      category=Category,
      color=get_color(category),
      module=LowGrip_Module,
      network_member=TRUE,
      x=runif(n(),-6,6),
      y=runif(n(),-40,40),
      z=sample(c(-1,1),n(),replace=TRUE)*runif(n(),18,28),
      size=.65,
      group="cross_associated"
    ) %>%
    select(
      id,label,category,color,module,network_member,
      x,y,z,size,group
    )
}else tibble()

# 其他283-node network疾病：更远的前后Z层
background_meta <- make_meta(background_codes)
n_background <- nrow(background_meta)

background_nodes <- if(n_background>0){
  background_meta %>%
    mutate(
      id=Disease_Code,
      label=Disease_Name,
      category=Category,
      color=get_color(category),
      module=LowGrip_Module,
      network_member=TRUE,
      x=runif(n(),-10,10),
      y=runif(n(),-55,55),
      z=sample(c(-1,1),n(),replace=TRUE)*runif(n(),35,60),
      size=.40,
      group="background_network"
    ) %>%
    select(
      id,label,category,color,module,network_member,
      x,y,z,size,group
    )
}else tibble()

# ── 12. 合并全部节点 ───────────────────────────────────────
nodes_3d <- bind_rows(
  center_node,
  up_pos,
  down_pos,
  up_assoc,
  down_assoc,
  cross_nodes_3d,
  background_nodes
)

duplicate_nodes <- nodes_3d %>%
  count(id,name="N") %>%
  filter(N>1)

if(nrow(duplicate_nodes)>0){
  print(duplicate_nodes)
  stop("❌ 节点ID重复")
}

nodes_3d <- nodes_3d %>%
  mutate(
    across(c(x,y,z,size),~round(.x,4))
  )

# ── 13. 完整LowGrip network edges ─────────────────────────
raw_network_edges <- igraph::as_data_frame(
  lowgrip_net,
  what="edges"
)

if("weight" %in% names(raw_network_edges)){
  background_network_edges <- raw_network_edges %>%
    transmute(
      from=clean_code(from),
      to=clean_code(to),
      weight=abs(as.numeric(weight)),
      type="network_background"
    )
}else{
  background_network_edges <- raw_network_edges %>%
    transmute(
      from=clean_code(from),
      to=clean_code(to),
      weight=.1,
      type="network_background"
    )
}

background_network_edges <- background_network_edges %>%
  filter(
    from %in% community_codes,
    to %in% community_codes,
    from!=to
  ) %>%
  mutate(pair=pair_key(from,to)) %>%
  distinct(pair,.keep_all=TRUE)

# ── 14. 原四层differential edges ──────────────────────────
association_up_edges <- fourth_up$edges %>%
  transmute(
    from=from,
    to=to,
    weight=.30,
    type="association_up"
  )

association_down_edges <- fourth_down$edges %>%
  transmute(
    from=from,
    to=to,
    weight=.30,
    type="association_down"
  )

# cross与up/down两侧真实differential边恢复
association_cross_edges <- bind_rows(
  fourth_up_all$edges %>%
    filter(to %in% cross_codes) %>%
    transmute(
      from=from,
      to=to,
      weight=.25,
      type="association_cross"
    ),
  fourth_down_all$edges %>%
    filter(to %in% cross_codes) %>%
    transmute(
      from=from,
      to=to,
      weight=.25,
      type="association_cross"
    )
) %>%
  distinct(from,to,.keep_all=TRUE)

# ── 15. 避免背景edge与高亮edge重复 ────────────────────────
highlight_pairs <- unique(c(
  pair_key(
    association_up_edges$from,
    association_up_edges$to
  ),
  pair_key(
    association_down_edges$from,
    association_down_edges$to
  ),
  pair_key(
    association_cross_edges$from,
    association_cross_edges$to
  )
))

background_network_edges <- background_network_edges %>%
  filter(!pair %in% highlight_pairs) %>%
  select(-pair)

# ── 16. Trajectory edges ──────────────────────────────────
trajectory_edges <- bind_rows(
  up_traj %>%
    transmute(
      from=Disease_Code,
      to=center_id,
      weight=abs(estimate),
      type="trajectory_up"
    ),
  down_traj %>%
    transmute(
      from=center_id,
      to=Disease_Code,
      weight=abs(estimate),
      type="trajectory_down"
    )
)

# ── 17. 合并全部edge ───────────────────────────────────────
edges_3d <- bind_rows(
  trajectory_edges,
  association_up_edges,
  association_down_edges,
  association_cross_edges,
  background_network_edges
) %>%
  distinct()

# ── 18. 最终严格核查 ───────────────────────────────────────
actual_disease_codes <- setdiff(
  unique(nodes_3d$id),
  center_id
)

missing_diseases <- setdiff(
  visual_disease_codes,
  actual_disease_codes
)

extra_diseases <- setdiff(
  actual_disease_codes,
  visual_disease_codes
)

missing_edge_nodes <- setdiff(
  unique(c(edges_3d$from,edges_3d$to)),
  nodes_3d$id
)

missing_community <- setdiff(
  community_codes,
  actual_disease_codes
)

cat("\n══════════════════════════════════════\n")
cat("             FINAL AUDIT\n")
cat("══════════════════════════════════════\n")
cat("Community diseases:",length(community_codes),"\n")
cat("Trajectory-only single-direction:",length(trajectory_only),"\n")
cat("Theoretical disease universe:",length(visual_disease_codes),"\n")
cat("Actual disease nodes:",length(actual_disease_codes),"\n")
cat("Total nodes including LowGrip:",nrow(nodes_3d),"\n")
cat("Missing diseases:",length(missing_diseases),"\n")
cat("Extra diseases:",length(extra_diseases),"\n")
cat("Missing Community diseases:",length(missing_community),"\n")
cat("Missing edge endpoints:",length(missing_edge_nodes),"\n")

cat("\nNode groups:\n")
print(nodes_3d %>% count(group))

cat("\nEdge groups:\n")
print(edges_3d %>% count(type))

if(length(missing_diseases)>0){
  cat("\n❌ Missing diseases:\n")
  print(missing_diseases)
}

if(length(extra_diseases)>0){
  cat("\n❌ Extra diseases:\n")
  print(extra_diseases)
}

if(length(missing_edge_nodes)>0){
  cat("\n❌ Missing edge endpoints:\n")
  print(missing_edge_nodes)
}

if(length(visual_disease_codes)!=286)
  stop("❌ 理论疾病universe不是286")

if(length(actual_disease_codes)!=286)
  stop("❌ 实际疾病节点不是286")

if(nrow(nodes_3d)!=287)
  stop("❌ 加LowGrip后总节点不是287")

if(length(missing_diseases)>0||
   length(extra_diseases)>0||
   length(missing_community)>0||
   length(missing_edge_nodes)>0)
  stop("❌ 完整性核查失败")

cat("\n✅ 所有完整性核查通过\n")

# ── 19. Audit文件 ─────────────────────────────────────────
node_audit <- nodes_3d %>%
  filter(id!=center_id) %>%
  transmute(
    Disease_Code=id,
    Disease_Name=label,
    Category=category,
    LowGrip_Module=module,
    Network_Member=network_member,
    Node_Group=group,
    x,y,z
  ) %>%
  arrange(
    factor(
      Node_Group,
      levels=c(
        "up_positive",
        "down_positive",
        "up_associated",
        "down_associated",
        "cross_associated",
        "background_network"
      )
    ),
    Disease_Code
  )

trajectory_only_audit <- make_meta(trajectory_only) %>%
  mutate(
    Direction=case_when(
      Disease_Code %in% up_codes~"Upstream",
      Disease_Code %in% down_codes~"Downstream",
      TRUE~"Unknown"
    ),
    Network_Member=FALSE
  )

bidirectional_audit <- make_meta(bidirectional) %>%
  mutate(
    In_Community=Disease_Code %in% community_codes,
    Display_Role=ifelse(
      In_Community,
      "Background network context",
      "Excluded from 286-disease visualization"
    )
  )

write_csv(
  node_audit,
  paste0(output_path,"Node_Audit_286Diseases.csv")
)

write_csv(
  trajectory_only_audit,
  paste0(output_path,"Trajectory_Only_Diseases.csv")
)

write_csv(
  make_meta(cross_codes),
  paste0(output_path,"Cross_Associated_Diseases.csv")
)

write_csv(
  make_meta(background_codes),
  paste0(output_path,"Background_Community_Diseases.csv")
)

write_csv(
  bidirectional_audit,
  paste0(output_path,"Bidirectional_Audit.csv")
)

write_csv(
  edges_3d,
  paste0(output_path,"Edge_Audit.csv")
)

# ── 20. JSON输出 ──────────────────────────────────────────
write_json(
  nodes_3d,
  paste0(output_path,"nodes.json"),
  pretty=TRUE,
  auto_unbox=TRUE,
  na="null"
)

write_json(
  edges_3d,
  paste0(output_path,"edges.json"),
  pretty=TRUE,
  auto_unbox=TRUE,
  na="null"
)

cat("\n══════════════════════════════════════\n")
cat("✅ JSON生成完成\n")
cat("Disease nodes:",length(actual_disease_codes),"\n")
cat("Total nodes:",nrow(nodes_3d),"\n")
cat("Total edges:",nrow(edges_3d),"\n")
cat("Output:",output_path,"\n")
cat("══════════════════════════════════════\n")