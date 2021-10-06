# Feedback in GUs
This git contains all the scripts and data used to build GENSOR Units from a regulatory network and calculate their feedback.

## Dependencies
- [pathway tools](https://biocyc.org/download.shtml)
- python3
- perl
- [perlcyc](https://solgenomics.net/downloads/perlcyc.pl) - perl module to remotely access pathway tools

## Usage 
Clone this repository locally
```
git clone https://github.com/georgette-femerling/Feedback-in-GUs.git
cd Feedback-in-GUs
```

### General use:
```
./Gensor-Units command [options]
```

### Updating GU library files
```
./Gensor-Units Update
```

### Building GUs from a Network
```
./Gensor-Units Gensor -n [path to network] -o [path to output directory]
```

### Finding Feedback
```
./Gensor-Units Feedback -g [path to Gensor Unit directory] -o [path to output directory]
```

## Notes on GENSOR units
Original scripts to build Gensor Units can be found in: [This Repository](https://github.com/dledezma/gensor_units)

## Authors 
- Gensor Unit builing - Daniela Ledezma-Tejeida
- Feedback Calculation and updates - Georgette Femerling-Romero
