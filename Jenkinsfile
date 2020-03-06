pipeline {
  agent {
    docker {
      image 'wolfram-docker-11.3.0'
    }

  }
  stages {
    stage('Run tests') {
      steps {
        dir(path: 'SpinWeightedSpheroidalHarmonics') {
          git 'https://github.com/BlackHolePerturbationToolkit/SpinWeightedSpheroidalHarmonics.git'
        }

        dir(path: 'KerrGeodesics') {
          git 'https://github.com/BlackHolePerturbationToolkit/KerrGeodesics.git'
        }

        dir(path: 'Teukolsky') {
          checkout scm
          sh 'Tests/AllTests.wls'
          junit 'TestReport.xml'
        }

      }
    }
  }
  options {
    timeout(time: 10, unit: 'MINUTES')
    skipDefaultCheckout(true)
  }
}