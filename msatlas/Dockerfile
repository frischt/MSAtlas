FROM r-base:3.5.2
MAINTAINER Tobias Frisch "frisch@imada.sdu.dk"

EXPOSE 3838

## install R dependencies
RUN apt-get update && apt-get install -y gdebi-core libssl-dev libcurl4-openssl-dev libssh2-1-dev libxml2-dev

#install java
RUN apt-get update && \
    apt-get install -y openjdk-8-jdk && \
    apt-get install -y ant && \
    apt-get install ca-certificates-java && \
    apt-get clean && \
    update-ca-certificates -f;

ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
RUN export JAVA_HOME

# more R dependencies
RUN R -e "install.packages(c('devtools'))"

# download shiny
RUN wget https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-1.5.7.907-amd64.deb && \
    gdebi -n shiny-server-1.5.7.907-amd64.deb && \
    rm shiny-server-1.5.7.907-amd64.deb

COPY install_depends.R install_depends.R
RUN R -e "source('install_depends.R')"

COPY . /srv/shiny-server

COPY shiny-server.sh /usr/bin/shiny-server.sh

CMD ["/usr/bin/shiny-server.sh"]
