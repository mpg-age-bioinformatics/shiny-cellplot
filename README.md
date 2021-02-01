## Cellplot shiny app

This app has been developed to work inside a docker container (Dockerfile file available [here](https://github.com/mpg-age-bioinformatics/shiny)).

To use this app locally you need start the container with:
```bash
docker run --rm -p 3838:3838 --name cellplot mpgagebioinformatics/shiny-cellplot
```
Access the app on your browser over [http://localhost:3838/cellplot](http://localhost:3838/cellplot).

The container can be stopped and the container removed with:
```bash
docker stop cellplot && docker rm cellplot
``` 
Removing the image once you've stopped the container:
```bash
docker rmi cellplot
```
