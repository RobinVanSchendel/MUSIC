plot_event_dna <- function(event, window = 20) {
  parsed <- str_split(event, "\\|", simplify = TRUE)
  events_df <- tibble(
    type  = parsed[1],
    start = as.numeric(parsed[2]),
    end   = as.numeric(parsed[3]),
    info = parsed[4]
  )
  
  events_df <- events_df %>%
    mutate(start = ifelse(type == "TINS",0,start),
           end = ifelse(type == "TINS",0,end),
           info = ifelse(type == "TINS","TINS",info),
           )
  
  ##helper function for MH lookup
  pos_to_index <- function(pos, window) {
    # If pos is length 0 or NA, return NA
    if (length(pos) == 0 || all(is.na(pos))) {
      return(NA_integer_)
    }
    pos + window + 1
  }
    
  #check if homology is present
  events_df <- events_df %>% mutate(mh_length = case_when(
    type == "PQ_DELETION" ~ gsub("bp","",info),
    type == "DELINS" ~ "0", ## all DELINS ## TODO: fix that
    type == "TINS" ~ "TINS", ## all DELINS ## TODO: fix that
    type == "HDR" ~ "0", ## all DELINS ## TODO: fix that
    grepl("DEL",type) ~ gsub("bp","",info), ## all deletions
    TRUE ~ info
  )) %>%
    mutate(mh_length = as.numeric(mh_length)) %>%
    mutate(mh_start = ifelse(mh_length>0,start-mh_length,NA))
    
  sequence = toupper("gcgtccaaggtcgggcaggaagagggcctatttcccatgat")
  # ---- sequence to positions ----
  bases <- str_split(sequence, "", simplify = TRUE)[1, ]
  seq_df <- tibble(
    pos  = -window:window,
    base = bases
  )
  
  ##add homology sequence to events_df
  mh_df <- events_df %>%
    filter(mh_length > 0)
  
  if(nrow(mh_df)>0){
    mh_df <- mh_df %>%
    rowwise() %>%
    mutate(
      mh_bases = list(
        bases[pos_to_index(mh_start, window) : pos_to_index(mh_start + mh_length - 1, window)]
      ),
      mh_pos = list(mh_start:(mh_start + mh_length - 1))
    ) %>%
    ungroup() %>%
    select(type, mh_bases, mh_pos) %>%
    tidyr::unnest(c(mh_bases, mh_pos)) %>%
    rename(base = mh_bases, pos = mh_pos)
  }
  
  # ---- handle 1-bp insertions ----
  ins_df <- events_df %>%
    filter( (grepl("INSERTION", type) | grepl("TINS",type) ) & info != "") 
  
  if(nrow(ins_df) > 0){
    ins_df <- ins_df %>%  # all insertions with info
    rowwise() %>%
    mutate(
      # Split multi-base insertions into individual letters
      base_list = list(strsplit(info, split = "")[[1]]),
      pos_list = list(seq(from = start, by = 1, length.out = length(base_list)))
    ) %>%
    ungroup() %>%
    tidyr::unnest(c(base_list, pos_list)) %>%
    rename(
      base = base_list,
      pos  = pos_list
    )
  }
  
  
  
  nudge = 0.5
  p <- ggplot() +
    
    # DNA sequence track
    geom_tile(
      data = seq_df,
      aes(x = pos, y = -1),
      height = 2,
      width  = 1,
      fill   = "white",
      color  = "white"
    ) +
    geom_text(
      data = seq_df,
      aes(x = pos, y = -1, label = base),
      family = "mono",
      size   = 5,
      fontface = "bold",
    ) +
    ##plot event
    geom_rect(data = events_df, aes(xmin = -window-nudge, xmax = start-nudge , ymin = 0, ymax = 1),fill = "#DADADA") +
    geom_rect(data = events_df, aes(xmin = end-nudge, xmax = window+nudge, ymin = 0, ymax = 1),fill = "#DADADA") +
    scale_x_continuous(
      limits = c(-window-nudge, window+nudge)
    )
    #add microhomology
    if(nrow(mh_df) > 0){
      p <- p +
      # Microhomology sequence inside rectangle
      geom_text(
        data = mh_df,
        aes(
          x = pos,
          y = 0.5,
          label = base
        ),
        fontface = "bold",
        family = "mono",
        size = 5,
        color = "steelblue"
      ) 
    }
    ## add insertions
    if(nrow(ins_df) > 0){
      p <- p +
        geom_text(
          data = ins_df,
          aes(x = pos, y = 0.5, label = base),  # y > 1 to place above event rectangle
          family = "mono",
          size = 5,
          fontface = "bold",
          color = "firebrick"
        )
    }
  
    ##cut site
    p <- p + geom_vline(xintercept = -.5, color = "red") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    labs(x = "Position", title = event)
    
    p
}
