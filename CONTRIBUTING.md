## Contributing

To help you with development we have deployed an RStudio server in the main container. 

After starting the container as descrived in the `README.md` you can start the RStudio server:
```bash
docker exec -i -t shiny sudo rstudio-server start
```
You can now access the RStudio server on you browser over [http://localhost:8787](http://localhost:8787) (user/password: mpiage/bioinf). This app can then be found in `/srv/shiny-server/cellplot`. If you need to add R libraries make sure you do so in the `requirements.R` file and that those libraries are added to the `/srv/shiny-server/cellplot/libs` folder. If you need to add system libraries you can do so by entering the container with:
```bash
docker exec -i -t shiny /bin/bash
```
Make sure any installed libraries are added to the end of the [Dockerfile](https://github.com/mpg-age-bioinformatics/shiny/blob/master/Dockerfile) and that the respective changes are committed.

Please make sure that all files downloaded by the user contain the version tag eg.:
```R
    # specify the output file name
    filename = function(){
      paste0('Histogram.',gitversion(),'.pdf')
    }
```

Stopping the RStudio server:
```bash
docker exec -i -t shiny sudo rstudio-server stop
```
Submitting your changes to the repo's app:
```bash
cd ~/shiny/cellplot
git add .
git commit -m "<describe your changes here>"
git push origin HEAD:master
```
Adding a tag:
```
git tag -e -a <major>.<minor>.<patch> HEAD
git push origin --tags
```
