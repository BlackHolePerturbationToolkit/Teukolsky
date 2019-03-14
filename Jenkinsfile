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

        dir(path: 'Teukolsky') {
          checkout scm
          sh 'Tests/AllTests.wls'
          junit 'TestReport.xml'
        }

      }
    }
  }
  options {
    skipDefaultCheckout(true)
  }
}