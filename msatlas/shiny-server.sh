#!/bin/sh

# Make sure the directory for individual app logs exists
mkdir -p /var/log/shiny-server
chown shiny.shiny /var/log/shiny-server

# Make sure the directory for session information is available
mkdir -p /home/atlas
chown shiny.shiny /home/atlas

exec shiny-server 2>&1
#exec shiny-server
