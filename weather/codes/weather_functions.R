
eto.penmon = function(tmax, tmin, wind_speed, solar_radiation, atm_pr, doy, lat, sealvl, vpd ){
  
  
  tmean = (tmax + tmin) / 2
  
  delta1 = (4098*(0.6108* exp(  (17.27*tmean)/(tmean +237.3))    ))/ ((tmean + 237.3)^2)
  
  ws = wind_speed
  
  tt = (900 / (tmean+273))*ws
  
  rs = solar_radiation
  
  pc = 0.000665*atm_pr #slope of the saturation vapor curve
  
  pt = pc  /  (  delta1 + pc*(1 + 0.34 * ws)  )
  
  dt = delta1 / (  delta1 + pc*(1 + 0.34 * ws)  )
  
  ea = 0.6108 * exp((17.27*tmin)/(tmin+237.3))
  
  dr = 1 + 0.033*cos(((2*pi)/365)*doy)
  
  delta = 0.409*sin(  ((2*pi/365)*doy)  - 1.39  )
  
  latrad = (pi/180)*lat
  
  ws = acos(-tan(latrad)*delta)
  
  ra = ((24*(60))/pi) * 0.082 * dr * ((ws*sin(latrad)*sin(delta))+(cos(latrad)*cos(delta)*sin(ws)))
  
  rns = (1- 0.23)*rs
  
  rso = (0.75 + 2e-5 * sealvl)*ra
  
  rnl = 4.903e-9*( (((tmax+273.16)^4)+ ((tmin+273.16)^4))/2 ) * (  0.34 - 0.14 * sqrt(ea) ) * (1.35 * (rs/rso) - 0.35)
  
  rn = rns - rnl
  
  rng = 0.408 * rn #missing an x here
  
  etrad = dt * rng
  etwind = pt * tt * vpd
  
  eto = (etrad + etwind)
  
  return(eto)
}

get_mesonet = function(data){
  data %>% 
    mutate_at(vars(ANTHESIS:HARVEST), ~as.Date(., format = '%m/%d/%Y')) %>% 
    #Express date in a format to read MESONET
    mutate_at(vars(SOWING, HARVEST), ~str_replace_all(as.character(as.Date(.,format = '%Y/%m/%d')), "-", "" ) ) %>% 
    group_by_at(vars(LOCATION:HARVEST)) %>% 
    nest() %>% 
    #Contruct a URL to access MESONET
    mutate(URL = paste0(# Location
      URL_, "?stn=", LOCATION_MESONET, 
      # Interval
      INTERVAL, 
      # Start Date
      "&t_start=", SOWING, "000000", 
      # Harvest
      "&t_end=", HARVEST, "000000", 
      #Variables to retrieve
      VARIABLES)
    ) %>% 
    # Get Mesonet
    mutate(mesonet = list(read.csv(URL))) %>% 
    ungroup() %>% 
    dplyr::select(LOCATION, mesonet) %>% 
    # Format of weather variables
    mutate(mesonet = mesonet %>% map(~.x %>% 
                                       mutate_at(vars(everything(),-TIMESTAMP, -STATION), ~as.numeric(.))
                                     )
           ) 
}

theme = theme(axis.title.y = element_text(size = rel(2.5), angle = 90), 
              axis.title.x = element_text(size = rel(2.5)), 
              #axis.text.x = element_text(hjust = 0.5, angle = 90, size=rel(2.0)), 
              #axis.text.y = element_text(hjust = 1, size=rel(2.0)), 
              axis.ticks = element_line(size = 0.7),
              axis.ticks.length = unit(0.4, "cm"),
              panel.background = element_rect(fill = "white",
                                              colour = "black",
                                              size = 0.5, linetype = "solid"),
              strip.text.x = element_text(size = 20),
              axis.text.y = element_text(hjust = 1, size=rel(2)), 
              axis.text.x = element_text(hjust = 0.5, angle = 0,  size=rel(2)),
              strip.background =element_rect(fill="white"),
              panel.border = element_rect(colour = "black", fill=NA),
              legend.title=element_text(size=rel(2.0)),
              legend.text=element_text(size=rel(2.0)), 
              #legend.position = c(0.15, 0.80),
              legend.background = element_rect(fill = "white", color = "black"), 
              #legend.position = 'top', 
              #legend.direction = "horizontal", 
              #legend.box = "horizontal"
)

