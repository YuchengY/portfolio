#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#

## mini-project2
## authors: Yucheng Yang, Dennis Wei

library(shiny)
library(rsconnect)
library(devtools)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(readr)
library(lubridate)
library(tidyr)
library(car)

#load the file
jurors<-read.csv("APM_jurors.csv")
trials<-read.csv("APM_trials.csv")
voir<-read.csv("APM_voir_dire_answers.csv")

#jurors and trials data sets joined by trial id
trialsJuror <- inner_join(trials, jurors, by = c("id" = "trial__id"))

#make factor var
trialsJuror$accepted <- as.factor(trialsJuror$struck_by)


#jurors and voir data sets joined by juror id
jurorsVoir <- inner_join(jurors, voir, by = c("id" = "juror_id"))


#three data sets joined together by juror id 
full<-full_join(trialsJuror,voir, by=c("id.y"= "juror_id"))



#data wrangling for interactive plot1
fig_1 <- trialsJuror %>%
    separate(cause_number, c("year","trial_num"),"-", convert=TRUE) %>%
    mutate(juror_race = race) %>%
    filter(juror_race == "Black" | juror_race == "White", accepted != "Juror excused/absent", accepted != "Unknown") %>%
    mutate(accepted = recode(accepted,"c('Struck by the defense',
                           'Struck by the state',
                           'Struck for cause',
                           'Struck without notation') = 'Struck';
                           c('Juror chosen to serve on jury',
                           'Juror not struck',
                           'Juror chosen as alternate') = 'Accepted'")) %>%
    
    group_by(year, accepted, race) %>%
    summarize(total=n()) 

##data wrangling for interactive plot2

table2 <- full %>%
  separate(cause_number, c("year","trial_num"),"-", convert=TRUE) %>%
  mutate(juror_race = race) %>%
  filter(defendant_race == "Black" | defendant_race == "White",
         juror_race == "Black" | juror_race == "White",
         accepted != "Juror excused/absent", accepted != "Unknown") %>%
  mutate(accepted = recode(accepted,"c('Struck by the defense',
                           'Struck by the state',
                           'Struck for cause',
                           'Struck without notation') = 'Struck';
                           c('Juror chosen to serve on jury',
                           'Juror not struck',
                           'Juror chosen as alternate') = 'Accepted'"), 
         new_verdict = ifelse(verdict=="Guilty on at least one offense", TRUE, FALSE)) 


# Define UI for application 
  
ui <- fluidPage(
  
  # Application title
  titlePanel("Visualizing jurors strike rates and race implications in Mississippi"),
  
  #tab setting 
  tabsetPanel(
    tabPanel("Plot1", fluid = TRUE,
             sidebarLayout(
               #add interactive slider and check box just to plot1 sidebar area 
               sidebarPanel(sliderInput(inputId = "slider", 
                                        label = "choose a year", 
                                        min = 1992, 
                                        max = max(fig_1$year), 
                                        value = 1998,
                                        step = 1,
                                        animate = TRUE),
                            br(),
                            checkboxInput(inputId = "color_check",
                                          label = "Check to color by race",
                                          value = TRUE)
               ),
               mainPanel(
                 plotOutput("plot1", click='plot_click'), dataTableOutput("table_1")
               )
             )
    ),
    tabPanel("Plot2", fluid = TRUE,
             sidebarLayout(
               #add interactive varSelectInput just to plot2 sidebar area 
               sidebarPanel(varSelectInput(inputId = "group", 
                                           label = "select your variable to group the data",
                                           data = select(table2, c(8,43,67,75,79,100)), 
                                           selected = "defendant_race"),
                            htmlOutput("text")),
               mainPanel(
                 plotOutput("plot2", click = 'plot_click2'),dataTableOutput("table_2")
                 
                 
               )
             )
    )
  )
)
    


server <- function(input, output) {
    
  #make plot showing the number of eligible jurors and allow user to choose the year shown
  #check box allow use to look at the number grouped by race, default is show grouping
    output$plot1 <- renderPlot({
        
        if (input$color_check){    
            ggplot(data = filter(fig_1, year <= input$slider), aes(x=year, y=total)) +  facet_wrap(~accepted)+
                geom_line(aes(color=race)) +
                geom_point(aes(color=race))+
                labs(title = "Figure 1: Number of Eligible Jurors Against Year",
                     x = "Year",
                     y = "Count",
                     color = "Race")}
        else{
            ggplot(data = filter(fig_1, year <= input$slider), aes(x=year, y=total)) +  facet_wrap(~accepted)+
                geom_line() +
                geom_point()+
                labs(title = "Figure 1: Number of Eligible Jurors Against Year",
                     x = "Year",
                     y = "Count"
                )
        }
    })
    
    output$plot2 <- renderPlot({
      
     
      #data wrangling to summarise struck rate grouped by year, juror_race and the variable user chose
      #and use ggplot to make line chart showing the time progression of struck rates
            table2 %>% group_by(year, juror_race, across(input$group)) %>%
                mutate(total = n(), struck = (accepted == "Struck")) %>%
                summarize(strike_rate = mean(struck)) %>% 
            ggplot(aes(x=year, y=strike_rate)) +
            geom_line(aes(color=juror_race)) +
            geom_point()+
            labs(title = str_c("Strike Rates Against Year grouped by ", input$group),
                 y = "Strike Rate",
                 x = "Years",
                 color = "Juror Race") +
            facet_wrap(as.formula(paste("~", input$group)))
      
      
      
        
    })   
   
    #making the output table for the point click in plot1
    output$table_1 <- renderDataTable({
        nearPoints(fig_1, input$plot_click)
    
  })
  
    #use html text output to explain to the users what those variables means
    output$text <- renderUI({
      str<- ""
      str_1<-"Notes on variables you can choose from:"
      str1 <- "defendant_race: the race of the defendant;"
      str2 <-  "gender: the gender of the juror;"
      str3 <- " crime_victim: juror has been the victim of a crime;"
      str4 <- "fam_law_enforcement: the juror has friends or family who work or have worked in law enforcement;"
      str5 <- " def_race: juror admitted the race of the defendant would affect his or her decision;"
      str6 <- " know_def: juror has prior familiarity with defendant through either personal or professional channels"
      
      HTML(paste(str, str_1, str1, str2, str3,str4,str5,str6, sep = '<br/>'))
   
})
    #making the output table for the point click in plot2
    output$table_2 <- renderDataTable({
      
      choose<-table2 %>% group_by(year, juror_race, across(input$group)) %>%
        mutate(total = n(), struck = (accepted == "Struck")) %>%
        summarize(strike_rate = mean(struck))
      
      nearPoints(choose, input$plot_click2)
    })

}
# Run the application 
shinyApp(ui = ui, server = server ,options = list(height =800))
