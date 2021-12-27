library(dplyr)
library(tidyverse)
library(readr)
library(lubridate)
library(tidyr)
library(car)
library(highcharter)
library(shiny)
library(caret)
library(rpart)
library(rpart.plot) 
library(shinydashboard)
library(shinyWidgets)
library(rsconnect)
library(e1071)

#reading in data
full_movies <- read.csv("full_movies.csv")

#recoding vars
full_movies <- full_movies %>% 
    mutate(content_rating = as.factor(content_rating)) %>%
    mutate(tomatometer_status = as.factor(tomatometer_status)) %>%
    mutate(audience_status = as.factor(audience_status)) %>%
    mutate(audience_count = log(audience_count)) %>%
    mutate(world_dollars = log(world_dollars)) %>%
    mutate(domestic_dollars = log(domestic_dollars)) %>%
    mutate(overseas_dollars = log(overseas_dollars))

#making another version without logged variables for prediction modelling
full_movies2<-read.csv("full_movies.csv")
full_movies2 <- full_movies2 %>% 
  mutate(content_rating = as.factor(content_rating)) %>%
  mutate(tomatometer_status = as.factor(tomatometer_status)) %>%
  mutate(audience_status = as.factor(audience_status)) %>%
  mutate(audience_count = log(audience_count)) 
full_movies_tree <- full_movies2 %>%  
    na_if("") %>% drop_na() %>% mutate()

#reading in data
withIMDB <- read_csv("withIMDB.csv")

#recoding vars
withIMDB <- withIMDB %>% 
    mutate(tomatometer_status = as.factor(tomatometer_status))

#ensure reproducibility
set.seed(45677)

#making training and test sets for the cross-validation 
n <- nrow(full_movies_tree)
train_index <- sample(nrow(full_movies_tree), size=round(.8*n)) # 80% training
movie_train <- full_movies_tree %>%slice(train_index)
movie_test <- full_movies_tree %>%slice(-train_index) 

#train control object
train_control <- trainControl(
    method =
        "cv"
    , # cross-validation
    number = 10 # 10-folds
)

#cross_validation to tune the cp for the optimal decision tree
rpart_cv <- train(
    tomatometer_status ~ runtime + content_rating + world_dollars + domestic_dollars + overseas_dollars,
    data = movie_train, # training data
    method ="rpart", # classification method
    tuneGrid = data.frame(cp = seq(0.0, .05, by = 0.0025)), # complexity parameters
    trControl = train_control # validation method
)
#the decision tree using the optimal cp
movie_rpart <- rpart(
    tomatometer_status ~ runtime + content_rating + world_dollars + domestic_dollars + overseas_dollars, 
    data = movie_train, cp=0.01)

#UI portion
ui <- fluidPage(
  headerPanel("Movie Analysis Based on IMDB and Rotten Tomatoes"),
  helpText("Authors: Cindy Guo, Dennis Wei, Yucheng Yang"),
    #setting aesthetics
    setBackgroundColor(
        color = c("AliceBlue"),
        gradient = c("linear", "radial"),
        direction = c("bottom", "top", "right", "left"),
        shinydashboard = FALSE
    ),
    tabsetPanel(
      tabPanel("About This Website", fluid = T,
               helpText("We are interested in exploring how movie characteristics have changed over the years, 
              and how these characteristics interact with each other.
              We looked at movie characteristics like domestic/overseas/world box office income, 
              movie ratings from Rotten Tomatoes and IMDB, 
              the number of ratings a movie gets, 
              movie runtimes and etc.
              We are interested in questions like:
              How do the movie ratings affect the box office income over the years?
              How do the domestic box office income and overseas box office income change over time? 
              When are the peaks and lows of film industry over the years?
              How does the movie length change over time?
              How do the audience ratings change overall?",
              br(),br(),
             
              strong("Note:"), "due to several large-scale variables 
              resulting in plots that are hard to intrepret,
              we log transformed variables of domestic/overseas/world 
              box office income and audience count when making 
              the corresponding plots.",
              br(),br(),
              strong( "Variables: "),
              br(),br(),
              "Runtime: The total runtime of each movie in minutes.",
              br(),br(),
              "Tomatometer rating: The tomatometer rating of each movie as found on the Rotten 
              Tomatoes website. The tomatometer ratings are derived from a mixture of certified critics 
              as well as the average tomatometer score, which represents the percentage of reviews positively 
              rating a movie. Generally, movies that rated at least 3.5 stars or better are 
              considered to be positively rated shows. 
              The ratings are separated into three categories. “Certified Fresh” movies have a tomatometer 
              score of 75% or above, and at least 5 top critic reviews. “Fresh” indicates that 60% of 
              reviewers gave the movie a positive rating. “Rotten” indicates that less than 60% of reviewers 
              gave the movie a positive review.",
              br(),br(),
              "Tomatometer count: The total number of certified rotten tomato critics.",
              br(),br(),
              "Audience rating: The audience rating is the percentage of users who have rated the show 
              positively, as defined above.",
              br(),br(),
              "Audience count: The total number of users who rated the movies. This variable has been logged 
              to better fit the graphs.",
              br(),br(),
              "World dollars: The total box office revenue of the movie. ",
              br(),br(),
              "Domestic dollars: The total box office revenue from the United states and Canada. ",
              br(),br(),
              "Overseas dollars: The total box office revenue from movies outside of the domestic dollar regions. ",
              br(),br(),
              "Imdb ratings: The average rating of imdb users. ",
              br(),br(),
              "Content rating: The minimum maturity level suggested to watch the movie. ",
              br(),br(),
              "Audience status: Spilled (less than 60% of users gave a rating of at least 3.5) or 
              Upright (at least 60% of users gave a rating of at least 3.5)")),
        #tab 1
        tabPanel("Continuous Variables Over Time",
                 sidebarPanel(
                     #set reactive facet widget
                     varSelectInput(inputId ="variablesy", label = "Variable Y:", 
                                    full_movies[c(12,15:16,18:19,24,25,27)]),
                     sliderInput(inputId = "slider", 
                                 label = "Choose a year", 
                                 min = min(full_movies$year), 
                                 max = max(full_movies$year), 
                                 value = max(full_movies$year),
                                 step = 1)
                 ),
                 mainPanel(tabsetPanel(
                    #outputs
                     tabPanel("Plot", plotOutput("plotOverall")),
                     tabPanel("Summary", plotOutput("plot", click = 'plot_click1'), dataTableOutput("table_1")),
                     tabPanel("Analysis", textOutput("analysis"))
                 )
                 )
        ), #tab 2
        tabPanel("Continuous Variables With Groupings",
                 tabPanel("Individual Movie", fluid = TRUE,
                                 sidebarLayout(
                                     sidebarPanel(
                                         #set reactive facet widget
                                         varSelectInput(inputId ="varTab1", label = "Variable Y:", 
                                                        full_movies[c(12,15:16,18:19,24,25,27)]),
                                         varSelectInput(inputId = "facet", 
                                                        label = "Select your variable to facet the plot",
                                                        data = select(full_movies, c(5)), 
                                                        selected = "content_rating"),
                                         varSelectInput(inputId = "choose_color", 
                                                        label = "Select your variable to color points",
                                                        data = select(full_movies, c(14,17)), 
                                                        selected = "tomatometer_status"),
                                         checkboxInput(inputId = "color2",
                                                       label = "Check to color",
                                                       value = FALSE),
                                         checkboxInput(inputId = "facet2",
                                                       label = "Check to facet",
                                                       value = FALSE)
                                       
                                     ),
                                     #output
                                     mainPanel(
                                         plotOutput("plot1")
                                     )
                                 )
                            )
                        ),
          #tab 3
          tabPanel("IMDB Ratings vs Rotton Tomatoes Ratings", fluid = TRUE,
               sidebarLayout(
                 sidebarPanel(
                   #set reactive facet widget
                   varSelectInput(inputId = "var_imdb", 
                                  label = "Select Y variable:",
                                  withIMDB[c(15,18)]),
                   checkboxInput(inputId = "imdb_facet",
                                 label = "Check to facet by tomatometer status",
                                 value = FALSE)
                 ),
                 #output
                 mainPanel(
                   plotOutput("plot_imdb"),
                   htmlOutput("text")
                 )
               )
      ),
      #tab 4
        tabPanel("Decision Tree",
                   tabPanel("Decision Tree", fluid = TRUE,
                            sidebarLayout(
                              sidebarPanel(
                                helpText("Explanation: A decision tree is a way to make predictions", 
                                         " based on a set of predictors. In this tab, you can ",
                                         "enter your movie runtime, content rating",
                                         " categories, and world/domestic/overseas box",
                                         " office income you expect,",
                                         " and we'll tell you the predicted Rotten Tomatoes status."),
                                #set reactive facet widget
                                wellPanel(
                                  textInput('runtime', 
                                                    "Enter runtime here:",
                                                    "200"),
                                        
                                          selectInput("content_rating", label = "enter content_rating value here", 
                                                      choices = list("PG" = "PG", "PG-13" = "PG-13", "R" = "R", "G" = "G", "NR" = "NR"), 
                                                      selected = "PG"),
                                        numericInput('world_dollars', "enter world_dollars here","50000000"),
                                        numericInput('domestic_dollars', "enter domestic_dollars here","25000000"),
                                        numericInput('overseas_dollars', "enter overseas_dollars here","25000000"),
    
                                          actionButton("submit","Submit")
                                      )
                                  ),
                                  #output
                                  mainPanel(
                                      plotOutput("plot_tree"),
                                      tableOutput('tree_table')
               )
          )
        )
      )
  )
)


# Define server logic required to draw a histogram
server <- function(input, output, session) {
    output$plot1 <- renderPlot({ 
        if (input$color2 & input$facet2){  
            ggplot(full_movies, aes_string(x=full_movies$year, y = input$varTab1, color=input$choose_color )) + 
                geom_jitter() +    
                facet_wrap(as.formula(paste("~", input$facet)))+
                stat_smooth(se=F, method = "lm")
            
        }else if (input$color2==F & input$facet2==T) {
            ggplot(full_movies, aes_string(x=full_movies$year, y = input$varTab1)) + 
                geom_jitter() +    
                facet_wrap(as.formula(paste("~", input$facet)))+
                stat_smooth(se=F, method = "lm")
 
        } else if (input$color2 == T & input$facet2 == F) {
            ggplot(full_movies, aes_string(x=full_movies$year, y = input$varTab1, color=input$choose_color)) + 
                geom_jitter() +  
                stat_smooth(se=F, method = "lm")
        } else {
            ggplot(full_movies, aes_string(x=full_movies$year, y = input$varTab1)) + 
                geom_jitter() +  
                stat_smooth(se=F, method = "lm")
        }
    })
    output$plot1_year <- renderPlot({ 
        full_movies_facet_year <- full_movies %>% 
            drop_na(!!input$varTab_year) %>% 
            group_by(year) %>% 
            summarise(mean = mean(eval(as.symbol(input$varTab_year)))) %>% 
            arrange(desc(mean))
        if (input$color2_year & input$facet2_year){  
            ggplot(full_movies_facet_year, aes_string(x=full_movies_facet_year$year, 
                                                      y = mean, 
                                                      color=input$choose_color_year )) + 
                geom_jitter() +    
                facet_wrap(as.formula(paste("~", input$facet_year)))+
                stat_smooth(se=F, method = "lm")
            
        }else if (input$color2_year==F & input$facet2_year==T) {
            ggplot(full_movies_facet_year, aes_string(x=full_movies_facet_year$year, 
                                                      y = mean)) + 
                geom_jitter() +    
                facet_wrap(as.formula(paste("~", input$facet_year)))+
                stat_smooth(se=F, method = "lm")
        } else if (input$color2_year == T & input$facet2_year == F) {
            ggplot(full_movies_facet_year, aes_string(x=full_movies_facet_year$year, 
                                           y = mean, 
                                           color=input$choose_color_year)) + 
                geom_jitter() +  
                stat_smooth(se=F, method = "lm")
        } else {
            ggplot(full_movies_facet_year, aes_string(x=full_movies_facet_year$year, 
                                                      y = mean)) + 
                geom_jitter() +  
                stat_smooth(se=F, method = "lm")
        }
    })
    #year~mean(box_office), select domestic/world/overseas box office
    output$plot3 <- renderPlot({ 
        full_movies %>% 
            drop_na(!!input$dollars) %>% 
            group_by(year) %>% 
            summarise(mean = mean(eval(as.symbol(input$dollars)))) %>% 
            arrange(desc(mean)) %>% 
            ggplot(aes(year, mean))+geom_point()+
            geom_line(color="blue", se=FALSE)+
            labs(title = str_c("mean ",  input$dollars, " by year"), y = str_c("mean ", input$dollars))
    })
    #output tab for click1
    output$table_1 <- renderDataTable({
        choose<-full_movies %>% 
            drop_na(!!input$variablesy) %>% 
            group_by(year) %>% 
            summarise(mean = mean(eval(as.symbol(input$variablesy)))) 
        nearPoints(choose, input$plot_click1)
    })
    summaryvar <- reactive({input$variablesy})
    #summary table for mean rotten_totaties/imdb ratings by year
    output$summary_table <- renderTable({
        full_movies %>% 
            drop_na(!!summaryvar()) %>% 
            group_by(year) %>% 
            summarise(mean = mean(eval(as.symbol(summaryvar())))) %>% 
            as_tibble()
    })
    output$plot <- renderPlot({
        full_movies_plot<- full_movies %>% 
            drop_na(!!input$variablesy) %>% 
            group_by(year) %>% 
            summarise(mean = mean(eval(as.symbol(input$variablesy)))) %>% 
            arrange(desc(mean)) 
        ggplot(data = filter(full_movies_plot, year <=input$slider), aes(year, mean))+geom_point()+
            geom_line(color="blue", se=FALSE)+
            labs(title = str_c("mean ",  input$variablesy, " by year"), 
                 y = str_c("mean ", input$variablesy),
                 x = "year")
    })
    #make plot using the year slider 
    output$plotOverall <- renderPlot({
        dataOverall <- filter(full_movies, year <=input$slider)
        ggplot(data = dataOverall, aes_string(x = dataOverall$year, 
                                              y = input$variablesy)) + 
            geom_point() +
            labs(title = str_c(input$variablesy, " by year"), 
                 y = str_c(input$variablesy),
                 x = "year")
    })
    #add analysis text 
    output$analysis <- renderText({
        if (input$variablesy == "runtime"){
            "The average runtime of movies fluctuated regularly throughout the years before settling down to an average runtime of 115 - 130 minutes 
              around the mid 1970s. This is quite interesting as it feels like in recent years the large blockbuster movies are becoming increasingly 
              long, for instance Avengers: Endgame was over 180 minutes long, but the average runtime is consistently within the range of 115-130 minutes. 
              It is interesting to note that the periods with lower mean runtime coincide with times of global issues. For instance, the 1940s coincide 
              with world war II, while 2020 coincides with the coronavirus pandemic. "
          
        } else if (input$variablesy == "overseas_dollars"){
            "Looking at the mean overseas box office income, 
            we see peaks at 1939 and 1942 as well (corresponding to the trend we see in domestic box office income and world box office income), 
            followed by a big drop in 1959 (since we only have one movies released in the 1959 in the data set ) where it hit the lowest point. 
            Then the number started to increase rapidly until 1973 when we see another peak. 
            Following that, we have a big drop again in 1976. 
            Then, we see an increase in that number in 1977 almost to the level it has during the highest peak in 1973. 
            The number gradually stabled ever since and maintained a relatively high level over the years. "
        } else if (input$variablesy == "domestic_dollars"){
            "Looking at the mean domestic box income for movies, 
            the first peak came in the year 1939, then 1973, 
            after that the highest mean in the history of the film industry occurs in the year 1983, 
            which is also when we had the “Star Wars: Episode VI - Return of the Jedi” released. 
            In the 80s, there were also many other 'blockbuster' films released like the Star War sequels.
            Ever since then, the mean domestic box office income was fluctuating declining."
            
        } else if (input$variablesy == "tomatometer_rating"){
            "The tomatometer rating trend suggests that the average rating of films, by rotten tomato critics,  has steadily decreased throughout 
              the years before slightly increasing again in recent years. There also appears to be an equal spread of “high” and “low” rated movies
              in any given year. These trends could simply reflect changes in how critics evaluate films, and how criteria might have changed over the years."
          
        } else if (input$variablesy == "tomatometer_count"){
            "The tomatometer count has steadily increased throughout the years, and experienced a sharp increase in the 2000s, perhaps owing to a 
              general increase in internet usage allowing for more people to become verified critics. This could also be attributed to an increase
              in movies produced giving users more opportunities to vote on movies. "
          
        } else if (input$variablesy == "audience_rating"){
            "The audience rating has generally followed the same overall trend as the tomatometer rating trends, but the overall rating is higher 
            suggesting that audiences tend to enjoy movies more than the “critics” might. "
          
        } else if (input$variablesy == "audience_count"){
            "The audience count on the other hand appears to follow a different pattern than the tomatometer count as it appears to have decreased over time. 
              Reasons for this shift might be due to more users becoming tomatometer critics, or perhaps migrating to other online movie platforms to rate movies. "
          
        } else if (input$variablesy == "world_dollars"){
            "The mean world box office income took a dump at 1960 (since we don't have movie released in the 1960 in the data set) from the peaks it has in 1939 and 1942, 
            hit its lowest point, then started to increase rapidly until the 70s and has its peaks in 1973, 
            1977, 1983. After that, the number started to stable down and was fluctuating increasing. 
            In 2020, possibly due to the Covid-19, the number decreased from the little peak it had in 2017 to the same level it had around 2000 and 2006."
          
        } 
    })
    #making the decision tree plot
    output$plot_tree<-renderPlot({
        rpart.plot(rpart_cv$finalModel)
    })
    #make world_dollar equal to overseas dollars added with domestic_dollar, 
    observeEvent(input$world_dollars,{
        new_domestic <-  input$world_dollars - input$overseas_dollars
        updateNumericInput(session, "domestic", value = new_domestic)
    })
    observeEvent(input$domestic_dollars,{
        newC <- input$world_dollars - input$domestic_dollars
        updateNumericInput(session, "overseas_dollars", value = newC)
    })
    observeEvent(input$overseas_dollars,{
        newA <- input$domestic_dollars+ input$overseas_dollars
        updateNumericInput(session, "world_dollars", value = newA)
    })
  #turn the user input into a data frame and predict outcome
    
    
    newdata <- reactive({
      
      test <- data.frame(runtime=input$runtime,content_rating=input$content_rating, 
                         world_dollars = input$world_dollars, domestic_dollars = input$domestic_dollars,
                         overseas_dollars= input$overseas_dollars) %>% mutate(content_rating = as.factor(content_rating), 
                                                                              runtime= as.integer(runtime), world_dollars = 
                                                                                as.numeric(world_dollars), domestic_dollars = 
                                                                                as.numeric(domestic_dollars), 
                                                                              overseas_dollars = as.numeric(input$overseas_dollars ))
      return(test)
      
    })
    
    output$tree_table <- renderTable({
        
        if (input$submit > 0) {
        preds <- predict(movie_rpart, newdata = newdata(), type ="class")
        newdata() %>%
            mutate(prediction = preds) %>% as_tibble()
        }
    })
    #making plots
    output$plot_imdb <- renderPlot({
        if (input$imdb_facet) {
            ggplot(withIMDB, aes_string(x = withIMDB$imdb_rating, 
                                        y = input$var_imdb)) +
                geom_jitter(size = 0.4) + 
                geom_smooth(method = "lm", se=F) + 
                facet_wrap(~tomatometer_status) + 
                labs(title = str_c(input$var_imdb, " against IMDB Ratings"), 
                     y = str_c(input$var_imdb),
                     x = "IMDB Ratings")
        } else {
            ggplot(withIMDB, aes_string(x = withIMDB$imdb_rating, 
                                        y = input$var_imdb)) +
                geom_jitter(size = 0.4) + 
                geom_smooth(method = "lm", se=F) +
                labs(title = str_c(input$var_imdb, " against IMDB Ratings"),
                     y = str_c(input$var_imdb),
                     x = "IMDB Ratings")
        }
    })
    #add text chunk for the plot
    output$text <- renderUI({
      str<- ""
      str1<- "Generally speaking, the imdb rating corresponds more positively to the audience rating"
      str2 <- " from rotten tomatoes. This can be explained by the fact that both sites are driven by"
      str3 <- "general audience members who aren’t necessarily approved “critics”. It is also worth noting" 
      str4 <-" that the imdb ratings and the rotten tomatoes ratings are more closely aligned in regards to"
      str5 <-" movies that are well-rated overall. On the other hand the rotten movies are negatively"
      str6 <- " correlated with the imdb ratings. However, this might be explained by a flaw in our dataset" 
      str7 <- " as our imdb ratings were gathered from the top 1000 rated imdb movies, thus all the imdb"
      str8 <- " movies only have positive ratings."
      HTML(paste(str, str1, str2, str3,str4,str5,str6, str7, str8, sep = '<br/>'))
    })
    
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
