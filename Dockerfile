FROM mpgagebioinformatics/shiny:0.4

RUN apt-get update && apt-get install -y git && apt-get clean && rm -rf /var/lib/apt/lists/*

ENV APP="cellplot"

RUN mkdir -p /srv/shiny-server/.git/modules/${APP}/refs/heads/

RUN git clone https://github.com/mpg-age-bioinformatics/shiny-${APP}.git /srv/shiny-server/${APP}

RUN cp /srv/shiny-server/${APP}/.git/refs/heads/master /srv/shiny-server/.git/modules/${APP}/refs/heads/