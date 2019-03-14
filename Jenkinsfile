pipeline {
  agent { docker { image 'wolfram-docker-11.3.0' } }
  options { skipDefaultCheckout true }
  stages {
    stage('Run tests') {
      steps {
        dir('SpinWeightedSpheroidalHarmonics') {
          git url: 'https://github.com/BlackHolePerturbationToolkit/SpinWeightedSpheroidalHarmonics.git'
        }
        dir('Teukolksy') {
          checkout scm
          sh 'Tests/AllTests.wls'
          junit 'TestReport.xml'
        }
      }
    }
  }
}