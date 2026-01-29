plot_event_dna <- function(event, window = 20) {
  parsed <- str_split(event, "\\|", simplify = TRUE)
  
  events_df <- tibble(
    type  = parsed[1],
    start = as.numeric(parsed[2]),
    end   = as.numeric(parsed[3]),
    info = parsed[4]
  )
  
  ##helper function for MH lookup
  pos_to_index <- function(pos, window) {
  pos + window + 1
  }
    
  #check if homology is present
  events_df <- events_df %>% mutate(mh_length = case_when(
    type == "PQ_DELETION" ~ gsub("bp","",info),
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
  
  nudge = 0.5
  ggplot() +
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
    ) +
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
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    labs(x = "Position", title = event)
}
