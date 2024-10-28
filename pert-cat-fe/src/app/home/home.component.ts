import { Component } from '@angular/core';
import {HeaderComponent} from "../shared/header/header.component";
import {MatCardModule} from "@angular/material/card";
import {FooterComponent} from "../shared/footer/footer.component";

@Component({
  selector: 'app-home',
  standalone: true,
  imports: [HeaderComponent, FooterComponent, MatCardModule],
  templateUrl: './home.component.html',
  styleUrl: './home.component.css'
})
export class HomeComponent {

}
